require(data.table)

args <- c("lag_effectiveness.rds","effectiveness.rds", "lag_comboeff.rds")
args <- commandArgs(trailingOnly = TRUE)
# args <- c("~/Dropbox/irs_timing-refit0_intro-fix.sqlite", "~/Dropbox/who/fig1_data/baseline.rds")

delay <- 2 # TODO could fetch this from db, and have on=.(year = year+delay)

combo.dt <- readRDS(args[1])
nolag_effectiveness.dt <- readRDS(args[2])

vecref <- nolag_effectiveness.dt[vc == 1 & vac == 0 & vc_coverage == 75]
vacref <- nolag_effectiveness.dt[vc == 0 & vac == 1 & vaccine == "traditional" & catchup == "catchup"]

veckeys <- setdiff(key(nolag_effectiveness.dt), c("vac","vc","vaccine","catchup"))
vackeys <- setdiff(key(nolag_effectiveness.dt), c("vac","vc","vc_coverage"))
bothkeys <- intersect(veckeys, vackeys)

vec.lag <- vecref[,
  .(vec.eff=eff, c.vec.eff=c.eff, year = year + delay),
  by=c(setdiff(veckeys,"year"))
][vacref, on=bothkeys, mult="first"][,
  .(vac_first = 1, vc_coverage=75,
    vec.eff=ifelse(is.na(vec.eff),0,vec.eff),
    c.vec.eff=ifelse(is.na(c.vec.eff),0,c.vec.eff),
    vac.eff = eff, c.vac.eff = c.eff,
    vaccine, catchup
    ),
  by=bothkeys
]

vac.lag <- vacref[,
  .(vac.eff=eff, c.vac.eff=c.eff, year = year + delay),
  by=c(setdiff(vackeys,"year"))
][vecref, on=bothkeys, mult="first"][,
.(vac_first = 0, vc_coverage,
 vac.eff=ifelse(is.na(vac.eff),0,vac.eff),
 c.vac.eff=ifelse(is.na(c.vac.eff),0,c.vac.eff),
 vec.eff = eff, c.vec.eff = c.eff,
 vaccine="traditional", catchup="catchup"
),
by=bothkeys
]

lags <- rbind(vec.lag, vac.lag)

syn.dt <- combo.dt[lags,
  on=.(vc_coverage, vaccine, catchup, vac_first, particle, replicate, year)
][, # get the interesting measures
  .(
    combo.eff=eff, ind.eff = (vec.eff + vac.eff - vec.eff*vac.eff),
    vec.eff, vac.eff,
    c.combo.eff=c.eff, c.ind.eff = (c.vec.eff + c.vac.eff - c.vec.eff*c.vac.eff),
    c.vec.eff, c.vac.eff
  ), keyby = c(union(vackeys, veckeys),"vac_first")
  # ...organized by relevant divisions
]

syn.dt[,
  syn := combo.eff - ind.eff
][,
  # if syn positive, should be fraction of space between 1 and ind.eff
  syn.frac := syn / ifelse(syn < 0, ind.eff, 1-ind.eff)
]

saveRDS(syn.dt, tail(args, 1))