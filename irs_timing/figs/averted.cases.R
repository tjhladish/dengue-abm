
args <- commandArgs(trailingOnly = TRUE)
# args <- c("~/Dropbox/who/fig1_data/stopping-baseline.cases.rds", "~/Dropbox/who/fig1_data/stopping-interventions.cases.rds", "~/Dropbox/who/fig1_data/stopping-averted.rds")
# args <- c("~/Dropbox/who/fig1_data/foi-baseline.cases.rds", "~/Dropbox/who/fig1_data/foi-interventions.cases.rds", "~/Dropbox/who/fig1_data/foi-averted.rds")

require(data.table)
# sessionInfo()

baseline.dt <- readRDS(args[1])
bkey <- grep("cases",names(baseline.dt),invert=T, value=T)
interventions.dt <- readRDS(args[2])
ikey <- grep("(cases|particle)",names(interventions.dt),invert=T, value=T)

avert <- interventions.dt[baseline.dt, on=bkey][,
  .(
    averted=i.cases-cases,
    cum.averted=i.cum.cases-cum.cases,
    i.cum.cases
  ),
  keyby=ikey
]

avert[, cum.eff := ifelse(i.cum.cases == 0 & cum.averted == 0, 0, cum.averted/i.cum.cases) ]

res <- avert[,
  .(
    averted.mn=mean(averted), averted.md=median(averted),
    cum.averted.mn=mean(cum.averted), cum.averted.md=median(cum.averted),
    cum.eff.mn=mean(cum.eff), cum.eff.md=median(cum.eff)
  ),
  by=ikey
]

saveRDS(res, args[3])
