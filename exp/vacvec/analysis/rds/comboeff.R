suppressPackageStartupMessages({
  require(data.table)
})

# debugging only lines; args overrides when actually making script
args <- c("effectiveness.rds", "comboeff.rds")
args <- c("foi_effectiveness.rds", "foi_comboeff.rds")
args <- commandArgs(trailingOnly = TRUE)

effectiveness.dt <- readRDS(args[1])

ekeys <- key(effectiveness.dt)

# extract the vector-control only results
veckeys <- setdiff(ekeys, c("scenario", "vaccine","catchup"))
vec.dt <- effectiveness.dt[scenario == "vc",
  .(vec.eff=eff, c.vec.eff=c.eff),
  keyby=veckeys
]

# extract the vaccine only results
vackeys <- setdiff(ekeys, c("scenario", "vc_coverage"))
vac.dt <- effectiveness.dt[scenario == "vac",
  .(vac.eff=eff, c.vac.eff=c.eff),
  keyby=vackeys
]
vac.dt[catchup=="vac-only", catchup := "vc+vac" ]

# extract the combo-only results
combo.dt <- effectiveness.dt[scenario == "vc+vac"]

syn.dt <- combo.dt[
  # join vector control only results
  vec.dt, on=veckeys, nomatch=0
][
  # join vaccine only results
  vac.dt, on=vackeys, nomatch=0
][,
  .(
    # joining to combo.dt, so eff / c.eff are the combo effs
    combo.eff = eff, c.combo.eff = c.eff,
    # as defined earlier when getting single interventions
    vec.eff, c.vec.eff, vac.eff, c.vac.eff,
    # based on naive independence assumption
    ind.eff = (vec.eff + vac.eff - vec.eff*vac.eff),
    c.ind.eff = (c.vec.eff + c.vac.eff - c.vec.eff*c.vac.eff)
  ), keyby = c(union(vackeys, veckeys))
]

syn.dt[,
  syn := combo.eff - ind.eff
][,
  # if syn positive, should be fraction of space between 1 and ind.eff
  syn.frac := syn / ifelse(syn < 0, ind.eff, 1-ind.eff)
]

saveRDS(syn.dt, tail(args, 1))