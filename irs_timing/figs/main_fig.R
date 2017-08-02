
require(ggplot2)
require(gridExtra)
require(data.table)
require(lubridate)

## get script args; for debugging, uncomment section that follows
args <- commandArgs(trailingOnly = TRUE)
# args <- c(
#   paste0("~/Dropbox/who/fig1_data/",
#          c("obs", "sim", "eip", "R0", "mos", "coverage", "duration", "durability"),
#          ".rds"),
#   "~/Dropbox/who/fig1_data/fig1.png"
# )

## read in assorted input data
obs.dt <- readRDS(args[1])
sim.dt <- readRDS(args[2])
eip.dt <- readRDS(args[3])
R0.dt  <- readRDS(args[4])
mos.dt <- readRDS(args[5])
coverage.dt <- readRDS(args[6])
duration.dt <- readRDS(args[7])
durability.dt <- readRDS(args[8])

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
doys <- 1:365-1 # need to make 0-364 to add a scale that starts w/ 1
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
  xmin=head(tick.locs[c(TRUE, FALSE)],-1), # every other month start, except last
  xmax=tick.locs[c(FALSE, TRUE)], # every other month start, except first
  ymin=-Inf, ymax=Inf, # setting these to Inf ensures the rects extend full width + don't mess w/ scale
  alpha=0.05 # adjust for how dark; ALT: might also change color?
)

# change this to change the base (plot) line thickness
# all other line thicknesses are set relative to this
ref.line.sz <- 1

baselinep <- ggplot(plot.dt, baseaes) +
  theme_minimal() +
  theme(
    legend.text.align = 0.0,
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(0.9, "cm"),
    legend.title = element_text(face="bold"),
    legend.direction = 'horizontal',
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
      `Mos. pop.`=expression(paste(M(t),'  ')),
      `R0`=expression(paste(R[0],'  ')),
      EIP="EIP",
      Effectiveness="Effectiveness",
      `Cases (observed)`="Cases (observed)",
      `Cases (model)`="Cases (model)   "
    )
  ) +
  scale_size_manual(
    values = c(reference=ref.line.sz,`150`=ref.line.sz*2,`90`=ref.line.sz,`30`=ref.line.sz/2),
    breaks = c("150","90","30"),
    name   = "IRS durability sensitivity",
    labels = function(x) sprintf("%3s day   ", x)
  ) +
  scale_x_continuous(
    # don't show the day-of-year scale - covered by the background month layer
    name = NULL, breaks = NULL, labels = NULL, limits = c(1,365), expand=c(0,0)
  ) +
  scale_linetype_manual(
    values = c(reference="solid",`1`="dashed",`90`="solid",`365`="dotted"),
    breaks = c("365","90","1"),
    name   = "IRS rollout sensitivity",
    labels = function(x) sprintf("%3s day   ", x)
  ) +
  guides(
    size = guide_legend(
      title.position = "top", title.hjust = 0
    ),
    color = guide_legend(
      override.aes = list(size=ref.line.sz), title.position = "top", title.hjust = 0
    ),
    linetype = guide_legend(
      override.aes = list(size=ref.line.sz), title.position = "top", title.hjust = 0
    )
  )

margin.theme <- function(t,r,b,l) theme(
  plot.margin = unit(c(t, r, b, l), "line")
)

small.gap <- 0.1
  big.gap <- small.gap * 2

# plot panel for month labels
p.month <- baselinep + annotate("text",
  x = name.locs, y = 0,
  label = c("J","F","M","A","M","J","J","A","S","O","N","D"), size = 10
) + scale_y_continuous(
  name=NULL, breaks=0, labels=NULL
) + theme(
  panel.grid = element_blank()
)

# legend for seasonal plots
seas.legend <- theme(
  legend.title = element_blank(),
  legend.position = c(0.15, 0.92), legend.justification = c(0, 0.9)
)

p.cases <- baselinep + geom_line(data=rbind(cases.dt)) +
  guides(size="none", linetype="none") +
  scale_y_continuous(name="Incidence", breaks = c(0,1), limits = c(0,1), labels = c(0,"Max")) +
  seas.legend

p.seasonal <- baselinep + geom_line(data=seasonal.dt[layer != "background" & duration == "reference"]) +
  geom_line(data=seasonal.dt[layer == "background" & duration == "reference"], size=0.7) +
  geom_hline(baseaes, seasonal.dt[duration == "365"]) +
  guides(size="none", linetype="none") +
  scale_y_continuous(name="Seasonal factors", breaks = c(0, 1), limits = c(0,1), labels = c(0,"Max")) +
  seas.legend

proactive.start <- yday(as_date("1970/6/1")) # June 1
proactive.end   <- proactive.start + 179 # campaign is 90 days, including day 1
reactive.start  <- yday(as_date("1970/11/1")) # Nov 1
reactive.end    <- yday(as_date("1970/11/1")+179)

pro.col <- "cyan4"
rea.col <- "darkgreen"

ln.size <- 5/3*ref.line.sz

# see https://github.com/tidyverse/ggplot2/pull/2132
# using absolute latest ggplot2 makes for sharp arrow
arrow_y_offset = 0.825
p.campaigns <- baselinep + annotate("segment",
  x = c(1,proactive.start)+1, xend=c(reactive.end,proactive.end), y = c(-arrow_y_offset,arrow_y_offset), yend=c(-arrow_y_offset,arrow_y_offset),
  color=c(rea.col,pro.col), size=ln.size, linejoin="mitre",#/5,
  arrow=arrow(35,unit(1.5,"line"),"last","closed")
) + annotate("segment",
  x = reactive.start+1, xend=365, y = -arrow_y_offset, yend = -arrow_y_offset,
  color=rea.col, size=ln.size
) + annotate("segment",
  x = c(reactive.start, proactive.start)+ln.size/2, xend=c(reactive.start, proactive.start)+ln.size/2, y = c(-1.5,1.5), yend=c(-0.15,0.15),
  color=c(rea.col, pro.col), size=ln.size
) + annotate("text",
  y=c(arrow_y_offset,-arrow_y_offset), x=c(proactive.start,reactive.start)-1,
  label=c("Proactive IRS", "Reactive IRS"), color=c(pro.col,rea.col), size = 10,
  hjust="right"
) + scale_y_continuous(
  name=NULL, breaks=0, limits = c(-1.5, 1.5), labels = NULL
) + theme(
  panel.grid = element_line(color = "white"),
  panel.background = element_rect(fill="grey85")
)

highlighter <- annotate("line",
  x=coverage.dt[coverage == 75 & duration == 90 & durability == 90 & layer == "foreground", doy],
  y=coverage.dt[coverage == 75 & duration == 90 & durability == 90 & layer == "foreground", value],
  size = ref.line.sz*3, color = "yellow"
)

legend.x <- 0.55

eff.legend <- theme(
  legend.position = c(legend.x,.94),
  legend.justification = c(0, 0.9)
)

#coverage.dt[,face:=""]

p.eff.coverage <- baselinep + #geom_line(data=coverage.dt) +
  highlighter + # facet_grid(face ~ .) +
  geom_line(aes(color=factor(coverage)), data=coverage.dt) +
  scale_color_manual(
    values = c(`25`="lightgrey",`50`="darkgrey",`75`="black"),
    breaks = c("75","50","25"),
    name = "IRS coverage sensitivity",
    labels = function(x) {
      c(sprintf("%s%%   ", head(x,-1)),sprintf("%s%%", tail(x,1)))
    }
  ) +
  guides(size="none", linetype="none") +
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(name="", limits = c(0,1)) +
  eff.legend

#duration.dt[,face:="Sensitivity analyses"]

p.eff.duration <- baselinep +
  highlighter + #facet_grid(face ~ .) +
  geom_line(data=duration.dt) +
  geom_hline(baseaes, duration.dt[duration == 365], show.legend = F) +
  guides(color="none", size="none") +
  scale_y_continuous(name="Effectiveness", limits = c(0,1), labels = function(x) sprintf("%1.2f", x)) +
  eff.legend

#durability.dt[,face:=""]

p.eff.durability <- baselinep +
  highlighter + #facet_grid(face ~ .) +
  geom_line(data=durability.dt) +
  guides(color="none", linetype="none") +
  scale_y_continuous(name="", limits = c(0,1)) +
  eff.legend

labeller <- function(l) annotate("text", x=10, y=.92, label=l, size=15)

mon.left  <- 4.88
mon.right <- 0.1 #1.85
seas.left <- .15
eff.right <- 0.1
eff.left  <- 0.1

# the final plotting arrangement
res = 300
mag = 0.85*res/72
png(args[9], width = 1000*mag, height = 1730*mag, units = "px", res=res)
grid.arrange(
  p.month                          + margin.theme(small.gap, mon.right, small.gap, mon.left),
  p.cases          + labeller("a") + margin.theme(small.gap, mon.right, small.gap,seas.left),
  p.seasonal       + labeller("b") + margin.theme(small.gap, mon.right, big.gap,seas.left),
  p.campaigns      + labeller("c") + margin.theme(big.gap,   mon.right, big.gap, mon.left),
  p.eff.coverage   + labeller("d") + margin.theme(big.gap,   eff.right, small.gap, eff.left),
  p.eff.duration   + labeller("e") + margin.theme(small.gap, eff.right, small.gap, eff.left),
  p.eff.durability + labeller("f") + margin.theme(small.gap, eff.right, small.gap, eff.left),
  p.month                          + margin.theme(small.gap, mon.right, small.gap,   mon.left),
  ncol=1,
  heights = c(0.1,0.6,0.6,0.2,1,1,1,0.1)
)
dev.off()
