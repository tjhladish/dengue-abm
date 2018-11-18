require(data.table)
require(ggplot2)

args <- c("baseline.rds","intervention.rds","effstats.rds","fig_2.png")
args <- commandArgs(trailingOnly = TRUE)

baseline.dt <- readRDS(args[1])
intervention.dt <- readRDS(args[2])
effstats.rds <- readRDS(args[3])

base.inc <- baseline.dt[,{
  qs <- quantile(s, probs = c(0.025, .25, .5, .75, .975))
  names(qs) <- c("lo.lo","lo","med","hi","hi.hi")
  as.list(qs)
}, keyby=.(vc_coverage = rep(0,length(year)), year)]

inte.inc <- intervention.dt[vac==0,{
  qs <- quantile(s, probs = c(0.025, .25, .5, .75, .975))
  names(qs) <- c("lo.lo","lo","med","hi","hi.hi")
  as.list(qs)  
}, keyby=.(vc_coverage, year)]

eff.mlt <- melt.data.table(
  effstats.rds[variable == "vec.eff", .(effectiveness=unique(med)), keyby=.(vc_coverage, year)],
  id.vars = c("vc_coverage","year"),
  variable.name = "measure"
)
  
plot.dt <- rbind(melt.data.table(
  rbind(base.inc, inte.inc)[,.(incidence=med),keyby=.(vc_coverage, year)],
  id.vars = c("vc_coverage","year"),
  variable.name = "measure"
), eff.mlt)

vec_lines <- c(`0`=3,`25`=2,`50`=5,`75`=1) # c(`25`="dotted",`50`="dashed",`75`="solid")

limits <- data.table(
  measure = factor(c(rep("incidence", 2), rep("effectiveness", 2))),
  vc_coverage = rep(25, 4),
  year = rep(1, 4),
  value = c(c(0., 22500.0), c(0.0, 1.0))
)

p<-ggplot(plot.dt) +
  aes(x=year, y=value, linetype=factor(vc_coverage)) +
  geom_line() +
  geom_blank(data=limits) +
  facet_grid(measure ~ ., scales="free", switch = "y") +
  scale_linetype_manual("VC Coverage %", values=vec_lines) +
  theme_minimal() + theme(
    legend.direction = "horizontal",
    axis.title.y = element_blank(),
    strip.placement = "outside",
    legend.position = c(0.5, 0.5),
    legend.justification = c(0.5, 0.5),
    legend.margin = margin(), legend.spacing = unit(0, "pt")
  )

ggsave(
  tail(args,1), p, device = "png",
  width = 6, height = 6, dpi = "retina", units = "in"
)