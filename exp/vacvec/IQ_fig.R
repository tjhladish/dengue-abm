require(data.table)
require(ggplot2)
# inter-quartile plots

args <- c("effectiveness.rds")
args <- commandArgs(trailingOnly = TRUE)

effectiveness.dt <- readRDS(args[1])

stats.dt <- effectiveness.dt[,{
  qs <- quantile(eff, probs = c(.25,.5,.75), na.rm = T)
  names(qs) <- c("IQlo","med","IQhi")
  as.list(qs)
}, keyby=.(vc, vac, vc_coverage, vac_mech, catchup, year)]

stats.dt[vc == 0, vc_coverage := 0]
stats.dt[, vac_mech := factor(ifelse(vac==1,c("cmdvi","trad")[vac_mech+1],"none"))]
stats.dt[, catchup := factor(c("none","catchup")[catchup+1])]

ggplot(stats.dt) + aes(x=year+1, y=med, ymin=IQlo, ymax=IQhi, fill=vac_mech, color=vac_mech) +
  facet_grid(catchup ~ vc_coverage) +
  geom_ribbon(mapping = aes(color=NULL), alpha=.5) + geom_line() +
  theme_minimal()