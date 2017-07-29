
require(ggplot2)
require(gridExtra)
require(data.table)
require(lubridate)

## get script args; for debugging, uncomment section that follows
args <- commandArgs(trailingOnly = TRUE)
# args <- c(
#   paste0("~/Dropbox/who/fig1_data/",
#          c("obs", "sim", "eip", "R0", "mos", "baseline", "interventions"),
#          ".rds"),
#   "~/Dropbox/who/fig1_data/fig1.png"
# )

## read in assorted input data
obs.dt <- readRDS(args[1])
sim.dt <- readRDS(args[2])
eip.dt <- readRDS(args[3])
R0.dt  <- readRDS(args[4])
mos.dt <- readRDS(args[5])
baseline.dt <- readRDS(args[6])
interventions.dt <- readRDS(args[7])

## perform effectiveness calcs
eff10.dt <- interventions.dt[baseline.dt, on="particle"][,
  # join baseline to interventions on particle basis
  # baseline has *only* particle as key
  .(eff10=(i.cases10-cases10)/i.cases10),
  keyby=.(doy, coverage, duration, durability, particle)
]

# take stats across particles; maintains separate results by intervention dimensions
stat.eff10.dt <- eff10.dt[,
  .(med.eff10 = stats::median(eff10)),
  keyby=.(doy, coverage, duration, durability)
]

running.mean.smooth.n <- 5
convo.factor <- rep(1, running.mean.smooth.n)/running.mean.smooth.n
# take doy-by-doy statistical results (med.eff10), and smooth with running mean (see ?filter)
stat.eff10.dt[duration != 365, # 365 duration has only one datum, no need (or ability) to smooth
  smooth := as.numeric(filter(
    med.eff10, convo.factor, method = "convolution",
    sides = 2, circular = T
    # want running mean to be centered (sides = 2) and wrap the series (circular = T)
  )),
  by=.(coverage, duration, durability)
]
stat.eff10.dt[duration == 365, smooth := med.eff10 ]

# convenience function to look at sensitivity studies by appropriate
# dimensions -- uses "filt" to fix some params + pickout correct other layers
# first item is raw data (value = med.eff10, layer = "background")
# second item is smoothed data (value = smooth, layer = "foreground")
slice <- function(filt) rbind(stat.eff10.dt[
  eval(filt), .(doy, value=med.eff10, variable="Effectiveness",
    coverage, duration, durability,
    layer = "background"
  )], stat.eff10.dt[
  eval(filt), .(doy, value=smooth, variable="Effectiveness",
    coverage, duration, durability,
    layer = "foreground"
)])

# set duration, durability, floating coverage
coverage.dt <- slice(expression(duration == 90 & durability == 90))
# set coverage, durability, floating duration
duration.dt <- slice(expression(coverage == 75 & durability == 90))
# set duration, coverage, floating durability
durability.dt <- slice(expression(coverage == 75 & duration == 90))

## assorted data.tables for plotting
cases.dt <- rbind(obs.dt, sim.dt) # panel a
seasonal.dt <- rbind(
  mos.dt[layer=="foreground"], # take out filter to show raw rainfall as well
  R0.dt,
  eip.dt
) # panel b
eff.dt <- rbind(coverage.dt, duration.dt, durability.dt) # panel c-e
plot.dt <- rbind(cases.dt, seasonal.dt, eff.dt) # combine all data to set universal scales

## build some reference doy-based locations for month columns, labels
doys <- 1:365-1
dates <- lubridate::as_date(doys) # convert doys into dates for use w/ lubridate; default year 1970 is non-leap
month.len <- rle(month.abb[month(dates)])$lengths # get the month lengths
tick.locs <- c(0, cumsum(month.len)) + 1 # place ticks at 0, plus the cumsum of each month length

# want: labels falling in middle of month boundaries
name.locs <- head(tick.locs,-1) + month.len/2 #

# overall approach to how we will show data dimensions
baseaes <- aes(
  x=doy, y=value, yintercept=value,
  # lines will day-of-year (1-365) against some value (10yr effectiveness, or scaled seasonal data)
  color=variable, # different colors for different kinds of things;
  # shape=factor(coverage), # exception: will use grey scale for coverage; tried shape, meh
  size=factor(durability),
  linetype=factor(duration),
  alpha=layer # used to show foreground (smoothed) vs background (raw) data
)

# use background highlight to indicate months
mon.cols <- annotate("rect",
  xmin=head(tick.locs[c(TRUE, FALSE)],-1),
  xmax=tick.locs[c(FALSE, TRUE)],
  ymin=-Inf, ymax=Inf,
  alpha=0.05
)

ref.line.sz <- 1

baselinep <- ggplot(plot.dt, baseaes) +
  theme_minimal() +
  theme(
    legend.text.align = 0.0,
    legend.key.width = unit(1, "cm"), legend.key.height = unit(0.9, "cm"),
    panel.grid.minor = element_blank(),
    text = element_text(size = 25),
    panel.background = element_rect(fill="grey98", color="white"),
    axis.title.y = element_text(margin=margin(r=10)),
    strip.text.y = element_text(angle=90)
  ) +
  mon.cols +
  scale_alpha_manual(values=c(foreground=1,background=0.3), guide="none") +
  scale_color_manual(values=c(
    `Mos. pop.`="blue",
    `R0`="#ff7700",
     `EIP`=3,
     `Effectiveness`="black",
     `Cases (observed)`="black",
     `Cases (model)`="red"
    ), labels=c(
      `Mos. pop.`="Mos. pop.",
      `R0`=expression(R[0]),
      EIP="EIP",
      Effectiveness="Effectiveness",
      `Cases (observed)`="Cases (observed)",
      `Cases (model)`="Cases (model)"
    )
  ) +
  scale_size_manual(
    values = c(reference=ref.line.sz,`150`=ref.line.sz*2,`90`=ref.line.sz,`30`=ref.line.sz/2),
    breaks = c("150","90","30"),
    name   = "IRS durability",
    labels = function(x) sprintf("%3s day", x)
  ) +
  scale_x_continuous(
    # don't show the day-of-year scale - covered by the background month layer
    name = NULL, breaks = NULL, labels = NULL, limits = c(1,365), expand=c(0,0)
  ) +
  scale_linetype_manual(
    values = c(reference="solid",`1`="dashed",`90`="solid",`365`="dotted"),
    breaks = c("365","90","1"),
    name   = "IRS rollout",
    labels = function(x) sprintf("%3s day", x)
  )

small.gap <- 0.1
  big.gap <- small.gap * 2

p.month <- baselinep + annotate("text",
  x = name.locs, y = 0,
  label = c("J","F","M","A","M","J","J","A","S","O","N","D"), size = 10
) + scale_y_continuous(
  name=NULL, breaks=0, labels=NULL
) + theme(
  panel.grid = element_blank()
)
mon.left <- 6.68
mon.right <- 2.25
p.month.top <- p.month + theme(plot.margin = unit(c(big.gap,   mon.right, small.gap, mon.left), "line"))
p.month.bot <- p.month + theme(plot.margin = unit(c(small.gap, mon.right, big.gap,   mon.left), "line"))

seas.legend <- theme(
  legend.title = element_blank(),
  legend.position = c(0.25, 0.92), legend.justification = c(0, 0.9)
)

ln.size.override <- guide_legend(override.aes = list(size=ref.line.sz))
line.override.col <- guides(color=ln.size.override)
line.override.lty <- guides(linetype=ln.size.override)

seas.left <- 1.95

p.cases <- baselinep + geom_line(data=rbind(cases.dt)) +
  guides(size="none", linetype="none") +
  scale_y_continuous(name="Incidence", breaks = c(0,1), limits = c(0,1), labels = c(0,"Max")) +
  seas.legend +
  theme(
    plot.margin = unit(c(small.gap,mon.right,small.gap,seas.left), "line")
  ) +
  line.override.col

p.seasonal <- baselinep + geom_line(data=seasonal.dt[layer != "background" & duration == "reference"]) +
  geom_line(data=seasonal.dt[layer == "background" & duration == "reference"], size=0.7) +
  geom_hline(baseaes, seasonal.dt[duration == "365"]) +
  guides(size="none", linetype="none") +
  scale_y_continuous(name="Seasonal factors", breaks = c(0, 1), limits = c(0,1), labels = c(0,"Max")) +
  seas.legend +
  theme(
    plot.margin = unit(c(small.gap,mon.right,big.gap,seas.left), "line")
  ) +
  line.override.col

proactive.start <- yday(as_date("1970/6/1")) # June 1
proactive.end   <- proactive.start + 179 # campaign is 90 days, including day 1
reactive.start  <- yday(as_date("1970/11/1")) # Nov 1
reactive.end    <- yday(as_date("1970/11/1")+179)

pro.col <- # "chocolate4"
rea.col <- "darkturquoise"

ln.size <- 5/3*ref.line.sz

p.campaigns <- baselinep + annotate("segment",
  x = c(1,proactive.start)+1, xend=c(reactive.end,proactive.end)-5, y = c(-.75,.75), yend=c(-.75,.75),
  color=c(rea.col,pro.col), size=ln.size#, arrow=arrow(35,unit(1.5,"line"),"last","closed")
) + annotate("segment",
          x = c(1,proactive.start)+1, xend=c(reactive.end,proactive.end), y = c(-.75,.75), yend=c(-.75,.75),
          color=c(rea.col,pro.col), size=ln.size/5,
          arrow=arrow(35,unit(1.5,"line"),"last","closed")
) + annotate("segment",
  x = reactive.start+1, xend=365, y = -.75, yend = -.75,
  color=rea.col, size=ln.size
) + annotate("segment",
  x = c(reactive.start, proactive.start)+ln.size/2, xend=c(reactive.start, proactive.start)+ln.size/2, y = c(-1.5,1.5), yend=c(0,0),
  color=c(rea.col, pro.col), size=ln.size
) + annotate("text",
  y=c(0.75,-0.75), x=c(proactive.start,reactive.start)-1,
  label=c("Proactive IRS", "Reactive IRS"), color='black', size = 10,
  hjust="right"
) + scale_y_continuous(name=NULL, breaks=0, limits = c(-1.5,1.5), labels = NULL) +
theme(
  panel.grid = element_line(color = "white"),
  panel.background = element_rect(fill="grey85"),
  plot.margin = unit(c(big.gap, mon.right, big.gap,   mon.left), "line")
)

legend.x <- 0.7

highlighter <- annotate("line",
  x=coverage.dt[coverage == 75 & duration == 90 & durability == 90 & layer == "foreground", doy],
  y=coverage.dt[coverage == 75 & duration == 90 & durability == 90 & layer == "foreground", value],
  size = ref.line.sz*3, color = "yellow"
)

eff.legend <- theme(
  legend.position = c(legend.x,.95),
  legend.justification = c(0,0.9),
  legend.title.align = 0.5,
  strip.background = element_rect(fill='black', color='black'),
  strip.text = element_text(color='white')
)

coverage.dt[,face:=""]

p.eff.coverage <- baselinep + #geom_line(data=coverage.dt) +
  highlighter + facet_grid(face ~ .) +
  geom_line(aes(color=factor(coverage)), data=coverage.dt) +
  scale_color_manual(
    values=c(`25`="lightgrey",`50`="darkgrey",`75`="black"),
    breaks=c("75","50","25"),
    name="IRS coverage",
    labels = function(x) sprintf("%s%%", x)
  ) +
  guides(size="none", linetype="none") +
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(name="", limits = c(0,1)) +
  eff.legend +
  theme(
    plot.margin = unit(c(big.gap,0.5,small.gap,1.9), "line")
  ) + line.override.col

duration.dt[,face:="Sensitivity analyses"]

p.eff.duration <- baselinep +
  highlighter + facet_grid(face ~ .) +
  geom_line(data=duration.dt) +
  geom_hline(baseaes, duration.dt[duration == 365], show.legend = F) +
  guides(color="none", size="none") +
  scale_y_continuous(name="Effectiveness", limits = c(0,1), labels = function(x) sprintf("%1.2f", x)) +
  eff.legend +
  theme(
    plot.margin = unit(c(small.gap,0.5,small.gap,1.9), "line")
  ) + line.override.lty

durability.dt[,face:=""]

p.eff.durability <- baselinep +
  highlighter + facet_grid(face ~ .) +
  geom_line(data=durability.dt) +
  guides(color="none", linetype="none") +
  scale_y_continuous(name="", limits = c(0,1)) +
  eff.legend +
  theme(
    plot.margin = unit(c(small.gap,0.5,small.gap,1.9), "line")
  )

labeller <- function(l) annotate("text", x=10, y=.92, label=l, size=15)

# the final plotting arrangement; TODO: move padding changes here to consolidate?
png(args[8], width = 1000, height = 1600, units = "px")
grid.arrange(
  p.month.top,
  p.cases          + labeller("a"),
  p.seasonal       + labeller("b"),
  p.campaigns      + labeller("c"),
  p.eff.coverage   + labeller("d"),
  p.eff.duration   + labeller("e"),
  p.eff.durability + labeller("f"),
  p.month.bot,
  ncol=1,
  heights = c(0.1,0.6,0.6,0.2,1,1,1,0.1)
)
dev.off()
