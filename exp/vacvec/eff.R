require(data.table)

args <- c("baseline.rds", "intervention.rds", "eff.rds")
args <- commandArgs(trailingOnly = TRUE)
# args <- c("~/Dropbox/irs_timing-refit0_intro-fix.sqlite", "~/Dropbox/who/fig1_data/baseline.rds")

baseline.dt <- readRDS(args[1])
intervention.dt <- readRDS(args[2])

## perform effectiveness calcs
eff.dt <- intervention.dt[baseline.dt, on=c("particle", "year")][,
  # join baseline to interventions on particle basis
  # baseline has *only* particle as key
  .(eff=(i.s-s)/i.s),
  keyby=.(vc, vac, vc_coverage, vac_mech, catchup, particle, year)
]

saveRDS(eff.dt, args[3])