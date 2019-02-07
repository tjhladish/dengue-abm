suppressPackageStartupMessages({
  require(data.table)
  require(cowplot)
})

args <- c("figref.rda", "rds/foi_v_mpop.rds", "fig/SIfig_7.png")
args <- commandArgs(trailingOnly = TRUE)

load(args[1])
bind.dt <- readRDS(args[2])

bind.dt[, measure := gsub("^.+([is])$","\\1", variable) ]
bind.dt[, context := gsub("^(.+)\\.[is]$","\\1", variable) ]
bind.dt[context %in% c("i","s"), context := "total"]

# stability check
# ggplot(bind.dt) + theme_minimal() + aes(x=year, y=med, alpha=foi, group=foi) + geom_line() + facet_grid(measure ~ context, scales="free_y")

plot.dt <- bind.dt[, .(y=mean(med)/yucpop, lo=mean(lo), hi=mean(hi)), keyby=.(foi, measure, context)]

p <- ggplot(plot.dt) +
  aes(x=foi, y=y, linetype=context) +
  facet_grid(measure ~ ., scales="free_y", labeller = labeller(
    measure=c(s="Cases",i="Infections")
  )) +
#  geom_ribbon(aes(y=NULL), alpha=0.5) +
  geom_line() +
  geom_limits(plot.dt[,.(y=c(0,max(y)), foi=0), by=measure]) +
  scale_y_continuous("Incidence per 100k", expand=c(0,0)) +
  scale_x_continuous("Relative Mosquito Population", expand=c(0,0)) +
  scale_linetype_manual("Source", labels=c(total="Both", intro="Introduced", local="Local"), values=c(total="solid", intro="dashed", local="dotted")) +
  coord_cartesian(clip = "off") +
  theme_minimal() + theme(
    panel.spacing = unit(18,"pt"),
    legend.position = c(.01,.99),
    legend.justification = c(0,1),
    strip.text.y = element_text(angle=90)
  )

# target rds should be last arg
# will sniff target name to determine proper keys


plotutil(p, h=6, w=6, tail(args, 1))