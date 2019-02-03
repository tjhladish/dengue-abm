suppressPackageStartupMessages({
  require(data.table)
  require(cowplot)
})

args <- c("figref.rda", "rds/foi_v_mpop.rds", "fig/SIfig_7.png")
args <- commandArgs(trailingOnly = TRUE)

load(args[1])
bind.dt <- readRDS(args[2])

# stability check
# ggplot(bind.dt) + aes(x=year, y=y, alpha=foi, group=foi) + geom_line() + facet_grid(measure ~ fraction)

p <- ggplot(bind.dt[, .(y=mean(y)), keyby=.(foi, measure, fraction)]) +
  aes(x=foi, y=log10(y), linetype=fraction) +
  facet_grid(measure ~ ., labeller = labeller(
    measure=c(cases="Cases",infections="All Infections")
  )) +
  geom_line() +
  scale_y_continuous("Log(Count)", expand=c(0,0)) +
  scale_x_continuous("Relative Mosquito Pop.", expand=c(0,0)) +
  scale_linetype_manual("Source", labels=c(total="Any", introduced="Introduced"), values=c(total="solid", introduced="dashed")) +
  coord_cartesian(xlim=c(0,1.5), ylim=c(2,6)) +
  theme_minimal() + theme(
    panel.spacing = unit(18,"pt")
  )

# target rds should be last arg
# will sniff target name to determine proper keys


plotutil(p, h=6, w=6, tail(args, 1))