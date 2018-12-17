require(data.table)

args <- c("baseline.rds", "intervention.rds", "eff.rds")
args <- commandArgs(trailingOnly = TRUE)
# args <- c("~/Dropbox/irs_timing-refit0_intro-fix.sqlite", "~/Dropbox/who/fig1_data/baseline.rds")

baseline.dt <- readRDS(args[1])
intervention.dt <- readRDS(args[2])

bkeys <- key(baseline.dt)
ikeys <- key(intervention.dt)

## perform effectiveness calcs
eff.dt <- intervention.dt[baseline.dt, on=bkeys, nomatch=0][,
  # join baseline to interventions on sample (e.g., particle; particle + replicate; p/r/foi; etc) basis
  .(
    bcases = i.s, icases = s,
    c.bcases = i.c.s, c.icases = c.s,
    eff = ifelse(i.s == s, 1.0, (i.s-s)/i.s),
    c.eff = ifelse(i.c.s == c.s, 1.0, (i.c.s-c.s)/i.c.s)
  ),
  keyby=ikeys
]

saveRDS(eff.dt, args[3])