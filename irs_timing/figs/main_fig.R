
require(ggplot2)
require(gridExtra)
require(data.table)
require(lubridate)

# args <- c(
#   paste0("~/Dropbox/who/fig1_data/",
#          c("obs", "sim", "eip", "R0", "mos", "baseline", "interventions"),
#          ".rds"),
#   "~/Dropbox/who/fig1_data/fig1.png"
# )
args <- commandArgs(trailingOnly = TRUE)
obs.dt <- readRDS(args[1])
sim.dt <- readRDS(args[2])
eip.dt <- readRDS(args[3])
R0.dt  <- readRDS(args[4])
mos.dt <- readRDS(args[5])
baseline.dt <- readRDS(args[6])
interventions.dt <- readRDS(args[7])

eff10.dt <- interventions.dt[baseline.dt, on="particle"][,
  .(eff10=(i.cases10-cases10)/i.cases10),
  keyby=.(doy, coverage, duration, durability, particle)
]

stat.eff10.dt <- eff10.dt[,
  .(med.eff10 = stats::median(eff10)),
  keyby=.(doy, coverage, duration, durability)
]

stat.eff10.dt[duration != 365,
  smooth := as.numeric(filter(
    med.eff10, rep(1/5,5), method = "convolution",
    sides = 2, circular = T
  )),
  by=.(coverage, duration, durability)
]
stat.eff10.dt[duration == 365, smooth := med.eff10 ]
# ggplot(stat.eff10.dt, aes(x=timing, y=med.eff10, linetype=factor(duration))) +
#   facet_grid(coverage ~ .) +
#   geom_line() + geom_hline(aes(yintercept=med.eff10, linetype=factor(duration)), stat.eff10.dt[duration == 2])

slice <- function(filt) rbind(stat.eff10.dt[
  eval(filt), .(doy, value=med.eff10, variable="Effectiveness",
    coverage, duration, durability,
    layer = "background"
  )], stat.eff10.dt[
  eval(filt), .(doy, value=smooth, variable="Effectiveness",
    coverage, duration, durability,
    layer = "foreground"
)])

coverage.dt <- slice(expression(duration == 90 & durability == 90))
duration.dt <- slice(expression(coverage == 75 & durability == 90))
durability.dt <- slice(expression(coverage == 75 & duration == 90))

cases.dt <- rbind(obs.dt, sim.dt)

seasonal.dt <- rbind(
  mos.dt[layer=="foreground"], # take out filter to show raw rainfall as well
  R0.dt,
  eip.dt
)

eff.dt <- rbind(coverage.dt, duration.dt, durability.dt)

plot.dt <- rbind(cases.dt, seasonal.dt, eff.dt)

doys <- 1:365-1
dates <- lubridate::as_date(doys) # convert doys into dates for use w/ lubridate; default year 1970 is non-leap
month.len <- rle(month.abb[month(dates)])$lengths # get the month lengths
tick.locs <- c(0, cumsum(month.len)) + 1 # place ticks at 0, plus the cumsum of each month length

# want: labels falling in middle of month boundaries
name.locs <- head(tick.locs,-1) + month.len/2 #

baseaes <- aes(
  x=doy, y=value, color=variable,
  size=factor(durability), alpha=layer, linetype=factor(duration),
  shape=factor(coverage), yintercept=value
)

baselinep <- ggplot(plot.dt, baseaes) +
  theme_minimal() +
  theme(
    legend.text.align = 0.0,
    legend.key.width = unit(1, "cm"), legend.key.height = unit(0.9, "cm"),
    panel.grid.minor = element_blank(),
    text = element_text(size = 25),
    plot.margin = unit(c(0.2,0.5,0.1,2),"line"),
    panel.background = element_rect(fill="grey98", color="white")
  ) +
  annotate("rect",
    xmin=head(tick.locs[c(TRUE, FALSE)],-1),
    xmax=tick.locs[c(FALSE, TRUE)],
    ymin=-Inf, ymax=Inf,
    alpha=0.05
  ) +
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
    values = c(reference=1,`150`=2,`90`=1,`30`=0.5),
    breaks = c("150","90","30"),
    name   = "IRS Durability Sensitivity",
    labels = function(x) sprintf("%3s day", x)
  ) +
  scale_x_continuous(
    name = NULL, breaks = NULL, labels = NULL, limits = c(1,365), expand=c(0,0)
  ) +
  scale_linetype_manual(
    values = c(reference="solid",`1`="dashed",`90`="solid",`365`="dotted"),
    breaks = c("365","90","1"),
    name   = "IRS Rollout Sensitivity",
    labels = function(x) sprintf("%3s day", x)
  )

p.month <- baselinep + annotate("text",
  x = name.locs, y = 0,
  label = c("J","F","M","A","M","J","J","A","S","O","N","D"), size = 10
) + scale_y_continuous(name="", breaks=0, labels="Mon") + theme(panel.grid = element_blank())

proactive.start <- yday(as_date("1970/6/1"))
proactive.end <- proactive.start + 179
reactive.start <- yday(as_date("1970/11/1"))
reactive.end <- yday(as_date("1970/11/1")+179)

p.campaigns <- baselinep + annotate("rect",
  xmin = c(1,proactive.start,reactive.start), xmax=c(reactive.end,proactive.end,365), ymin = c(-0.5,0,-0.5), ymax=c(0,0.5,0),
  fill=c("blue","green","blue"), alpha = 0.5
) + annotate("text",
  y=c(0.25,-0.25), x=c((proactive.start+proactive.end)/2,(1+reactive.end)/2),
  label=c("proactive IRS", "reactive IRS"), size = 10
) + scale_y_continuous(name=NULL, breaks=0, labels = "Timing") + theme(panel.grid = element_blank())

labeller <- function(l) annotate("text", x=15, y=1, label=l, size=15)

p.cases <- baselinep + geom_line(data=rbind(cases.dt)) +
  guides(size="none", linetype="none") +
  scale_y_continuous(name="Incidence", breaks = c(0,1), labels = c(0,"Max")) +
  theme(legend.title = element_blank(), legend.position = c(0.1,0.85), legend.justification = c(0,0.9)) +
  labeller("a")

p.seasonal <- baselinep + geom_line(data=seasonal.dt[layer != "background"]) +
  geom_line(data=seasonal.dt[layer == "background"], size=0.7) +
  guides(size="none", linetype="none") +
  scale_y_continuous(name="Seasonal Factors", breaks = c(0,1), labels = c(0,"Max")) +
  theme(legend.title = element_blank(), legend.position = c(0.2,0.85), legend.justification = c(0,0.9)) +
  labeller("b")

legend.x <- 0.7

highlighter <- annotate("line",
  x=coverage.dt[coverage == 75 & duration == 90 & durability == 90 & layer == "foreground", doy],
  y=coverage.dt[coverage == 75 & duration == 90 & durability == 90 & layer == "foreground", value],
  size = 3, color = "yellow", alpha = 0.5
)

eff.legend <- theme(
  legend.position = c(legend.x,.95),
  legend.justification = c(0,0.9),
  legend.title.align = 0.5
)

p.eff.coverage <- baselinep + #geom_line(data=coverage.dt) +
  highlighter +
  geom_line(aes(color=factor(coverage)), data=coverage.dt) +
  scale_color_manual(
    values=c(`25`="lightgrey",`50`="darkgrey",`75`="black"),
    breaks=c("75","50","25"),
    name="IRS Coverage Sensitivity",
    labels = function(x) sprintf("%s%%", x)
  ) +
  guides(size="none", linetype="none") +
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(name="", limits = c(0,1)) +
  eff.legend + labeller("c")

p.eff.duration <- baselinep +
  highlighter +
  geom_line(data=duration.dt) +
  geom_hline(baseaes, duration.dt[duration == 365], show.legend = F) +
  guides(color="none", size="none") +
  scale_y_continuous(name="Effectiveness", limits = c(0,1), labels = function(x) sprintf("%1.2f", x)) +
  eff.legend + labeller("d")

p.eff.durability <- baselinep +
  highlighter +
  geom_line(data=durability.dt) +
  guides(color="none", linetype="none") +
  scale_y_continuous(name="", limits = c(0,1)) +
  eff.legend + labeller("e")

png(args[8], width = 1000, height = 1600, units = "px")
grid.arrange(
  p.month,
  p.cases, p.seasonal,
  p.campaigns,
  p.eff.coverage, p.eff.duration, p.eff.durability,
  p.month,
  ncol=1,
  heights = c(0.1,0.6,0.6,0.2,1,1,1,0.1)
)
dev.off()
