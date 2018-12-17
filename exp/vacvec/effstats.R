require(data.table)

args <- c("comboeff.rds", "effstats.rds")
args <- c("foi_comboeff.rds", "foi_effstats.rds")
args <- c("lag_comboeff.rds", "lag_effstats.rds")
args <- commandArgs(trailingOnly = TRUE)
# args <- c("~/Dropbox/irs_timing-refit0_intro-fix.sqlite", "~/Dropbox/who/fig1_data/baseline.rds")

comboeff.dt <- readRDS(args[1])
ckeys <- key(comboeff.dt)

mlt <- melt.data.table(
  comboeff.dt,
  id.vars = ckeys
)

mkeys <- c("variable", setdiff(ckeys,c("particle","replicate")))

res <- mlt[,{
  qs <- quantile(value, probs = c(0.025,.25,.5,.75,.975), na.rm=T)
  names(qs) <- c("lo.lo","lo","med","hi","hi.hi")
  as.list(qs)
  # .(med=qs[3])
}, keyby=mkeys]

saveRDS(res, args[2])