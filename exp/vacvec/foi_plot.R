require(data.table)
require(ggplot2)

args <- c("foi_baseline.rds","foi_intervention.rds","foi_effectiveness.rds", "foi_plot.png")
args <- commandArgs(trailingOnly = TRUE)

baseline.dt <- readRDS(args[1])
intervention.dt <- readRDS(args[2])
effstats.rds <- readRDS(args[3])

bkeys <- setdiff(key(baseline.dt), c("particle","replicate"))

base.inc <- baseline.dt[,{
  qs <- quantile(s, probs = c(0.025, .25, .5, .75, .975))
  names(qs) <- c("lo.lo","lo","med","hi","hi.hi")
  c(as.list(qs), vc_coverage=0)
}, keyby=bkeys]

scenkeys <- setdiff(key(intervention.dt), c("particle","replicate"))

inte.inc <- intervention.dt[,{
  qs <- quantile(s, probs = c(0.025, .25, .5, .75, .975))
  names(qs) <- c("lo.lo","lo","med","hi","hi.hi")
  as.list(qs)  
}, keyby=scenkeys]

yucpop <- 18.17734

ggplot(inte.inc) + aes(x=year, y=med/yucpop, linetype=vaccine, size=factor(vc_coverage)) +
  facet_grid(foi ~ .) + geom_line() +
  scale_size_manual(values=c(`0`=.75, `75`=1.5))

eff.mlt <- melt.data.table(
  effstats.rds[variable == "vec.eff", .(Effectiveness=unique(med)), keyby=.(vc_coverage, year)],
  id.vars = c("vc_coverage","year"),
  variable.name = "measure"
)
  
plot.dt <- rbind(melt.data.table(
  rbind(base.inc, inte.inc)[,.(Incidence=med),keyby=.(vc_coverage, year)],
  id.vars = c("vc_coverage","year"),
  variable.name = "measure"
), eff.mlt)

vec_lines <- c(`0`=3,`25`=2,`50`=5,`75`=1) # c(`25`="dotted",`50`="dashed",`75`="solid")

limits <- data.table(
  measure = factor(c(rep("Incidence", 2), rep("Effectiveness", 2))),
  vc_coverage = rep(25, 4),
  year = rep(1, 4),
  value = c(c(0., 22500.0), c(0.0, 1.0))
)

p<-ggplot(plot.dt) +
  aes(x=year, y=value, linetype=factor(vc_coverage)) +
  geom_line() +
  geom_blank(data=limits) +
  facet_grid(measure ~ ., scales="free", switch = "y") +
  scale_linetype_manual("Coverage %", values=vec_lines) +
  scale_x_continuous("Year", expand = expand_scale(0.03, 0)) +
  scale_y_continuous(expand = expand_scale(0, 0)) +
  theme_minimal() + theme(
    legend.direction = "horizontal",
    axis.title.y = element_blank(),
    strip.placement = "outside",
    legend.position = c(0.5, 0.5),
    legend.justification = c(0.5, 0.5),
    legend.margin = margin(), legend.spacing = unit(0, "pt"),
    panel.spacing.y = unit(15, "pt")
  )

ggsave(
  tail(args,1), p, device = "png",
  width = 4, height = 6, dpi = "retina", units = "in"
)