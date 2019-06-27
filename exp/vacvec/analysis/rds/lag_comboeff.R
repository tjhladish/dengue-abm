suppressPackageStartupMessages(
	require(data.table)
)

args <- c("all_effectiveness.rds", "lag_comboeff.rds")
args <- commandArgs(trailingOnly = TRUE)
# args <- c("~/Dropbox/irs_timing-refit0_intro-fix.sqlite", "~/Dropbox/who/fig1_data/baseline.rds")

combo.dt <- readRDS(args[1])[ivn_lag != 0]
nolag_effectiveness.dt <- readRDS(args[1])[ivn_lag == 0 & foi == 1.0]
allkeys <- key(nolag_effectiveness.dt)

delays <- combo.dt[,unique(ivn_lag)]

delayer <- function(del, dt) { 
  slice <- dt[order(year),
    .(year, bcases, eff=c(rep(0, del), head(eff,-del)), c.bcases, ivn_lag = del),
    keyby=.(particle, replicate)
  ]
  slice[, icases := ceiling((1-eff)*bcases) ][
    order(year), c.icases := cumsum(icases), by=.(particle, replicate)
  ]
  slice[, c.eff := ifelse(c.bcases == c.icases, 0, (c.bcases - c.icases)/c.bcases) ]
  return(slice)
}

vecnolag <- nolag_effectiveness.dt[scenario == "vc" & vc_coverage == 75]
vacnolag <- nolag_effectiveness.dt[scenario == "vac" & vaccine == "d70e" & catchup == "vac-only"]

veclag <- rbindlist(lapply(delays, delayer, dt=vecnolag))[, .(vec.eff=eff, c.vec.eff=c.eff),keyby=.(ivn_lag, particle, replicate, year)]
vaclag <- rbindlist(lapply(delays, delayer, dt=vacnolag))[, .(vac.eff=eff, c.vac.eff=c.eff),keyby=.(ivn_lag, particle, replicate, year)]

vacref <- vacnolag[,.(vac.eff = eff, c.vac.eff = c.eff), keyby=.(particle, replicate, year)]
vecref <- vecnolag[,.(vec.eff = eff, c.vec.eff = c.eff), keyby=.(particle, replicate, year)]

vec.lag <- veclag[vacref, on=key(vacref)][, vac_first := 1 ]
vac.lag <- vaclag[vecref, on=key(vecref)][, vac_first := 0 ]

lags <- rbind(vec.lag, vac.lag)

syn.dt <- combo.dt[lags,
  on=.(vac_first, ivn_lag, particle, replicate, year),
  nomatch=0
][, # get the interesting measures
  .(
    combo.eff=eff, ind.eff = (vec.eff + vac.eff - vec.eff*vac.eff),
    vec.eff, vac.eff,
    c.combo.eff=c.eff, c.ind.eff = (c.vec.eff + c.vac.eff - c.vec.eff*c.vac.eff),
    c.vec.eff, c.vac.eff
  ), keyby = c("vac_first", "ivn_lag", allkeys)
  # ...organized by relevant divisions
]

syn.dt[,
  syn := combo.eff - ind.eff
][,
  # if syn positive, should be fraction of space between 1 and ind.eff
  syn.frac := syn / ifelse(syn < 0, ind.eff, 1-ind.eff)
]

saveRDS(syn.dt, tail(args, 1))