require(data.table)

args <- c("comboeff.rds", "effstats.rds")
args <- commandArgs(trailingOnly = TRUE)
# args <- c("~/Dropbox/irs_timing-refit0_intro-fix.sqlite", "~/Dropbox/who/fig1_data/baseline.rds")

comboeff.dt <- readRDS(args[1])

mlt <- melt.data.table(
  comboeff.dt,
  id.vars = c("year","vc_coverage","vac_mech","catchup", "particle", "replicate")
)

setkeyv(mlt, c("variable", "vc_coverage","vac_mech", "catchup", "year"))

res <- mlt[,{
  qs <- quantile(value, probs = c(0.025,.25,.5,.75,.975), na.rm=T)
  names(qs) <- c("lo.lo","lo","med","hi","hi.hi")
  as.list(qs)
  # .(med=qs[3])
}, by=key(mlt)]

saveRDS(res, args[2])