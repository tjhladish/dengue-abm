
args <- commandArgs(trailingOnly = TRUE)
# args stopping-baseline.rds stopping-interventions.rds averted.cases.rds

require(data.table)

baseline.dt <- readRDS(args[1])
interventions.dt <- readRDS(args[2])

avert <- interventions.dt[baseline.dt, on=.(particle, year)][,
  .(
    averted=i.cases-cases,
    cum.averted=i.cum.cases-cum.cases,
    i.cum.cases
  ),
  keyby=.(coverage, particle, year)
]

avert[, cum.eff := ifelse(i.cum.cases == 0 & cum.averted == 0, 0, cum.averted/i.cum.cases) ]

res <- avert[,
  .(
    averted.mn=mean(averted), averted.md=median(averted),
    cum.averted.mn=mean(cum.averted), cum.averted.md=median(cum.averted),
    cum.eff.mn=mean(cum.eff), cum.eff.md=median(cum.eff)
  ),
  by=.(coverage, year)
]

saveRDS(res, args[3])