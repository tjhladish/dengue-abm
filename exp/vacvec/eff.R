require(data.table)

args <- c("baseline.rds", "intervention.rds", "eff.rds")
args <- commandArgs(trailingOnly = TRUE)
# args <- c("~/Dropbox/irs_timing-refit0_intro-fix.sqlite", "~/Dropbox/who/fig1_data/baseline.rds")

baseline.dt <- readRDS(args[1])
intervention.dt <- readRDS(args[2])

ikeys <- key(intervention.dt)

## perform effectiveness calcs
eff.dt <- intervention.dt[baseline.dt, on=.(particle = particle, replicate = replicate, year = year), nomatch=0][,
  # join baseline to interventions on particle basis
  # baseline has *only* particle as key
  .(
    eff = ifelse(i.s == s, 1.0, (i.s-s)/i.s),
    c.eff = ifelse(i.c.s == c.s, 1.0, (i.c.s-c.s)/i.c.s)
  ),
  keyby=ikeys
]

saveRDS(eff.dt, args[3])