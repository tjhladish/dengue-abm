require(data.table)

args <- c("effectiveness.rds", "comboeff.rds")
args <- commandArgs(trailingOnly = TRUE)
# args <- c("~/Dropbox/irs_timing-refit0_intro-fix.sqlite", "~/Dropbox/who/fig1_data/baseline.rds")

effectiveness.dt <- readRDS(args[1])

vec.dt <- effectiveness.dt[vac == 0, .(vec.eff=eff), by=.(particle, replicate, vc_coverage, year)]
vac.dt <- effectiveness.dt[vc == 0, .(vac.eff=eff), by=.(particle, replicate, vac_mech, catchup, year)]
combo.dt <- effectiveness.dt[(vc != 0) & (vac != 0)]

syn.dt <- combo.dt[
  vec.dt, on=.(particle, replicate, vc_coverage, year), nomatch=0 # join vector control only results
][
  vac.dt, on=.(particle, replicate, vac_mech, catchup, year), nomatch=0 # join vaccine only results
][, # get the interesting measures
  .(
    combo.eff=eff, ind.eff = (vec.eff + vac.eff - vec.eff*vac.eff),
    vec.eff, vac.eff
  ), by=.(particle, replicate, year, vc_coverage, vac_mech, catchup)
  # ...organized by relevant divisions
]

syn.dt[,
  syn := combo.eff - ind.eff
][,
  # if syn positive, should be fraction of space between 1 and ind.eff
  syn.frac := syn / ifelse(syn < 0, ind.eff, 1-ind.eff)
]

saveRDS(syn.dt, args[2])

stat.syn.dt <- syn.dt[,
  .(
    med.syn = stats::median(syn, na.rm=T),
    med.syn.frac = stats::median(syn.frac, na.rm=T)
  ),
  #.(med.eff = sum(eff)),
  keyby=.(
    vc_coverage,
    vac_mech,
    catchup,
    year
  )
]