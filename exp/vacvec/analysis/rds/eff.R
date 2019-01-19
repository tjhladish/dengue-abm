suppressPackageStartupMessages({
  require(data.table)
})

# development
args <- c("baseline.rds", "intervention.rds", "effectiveness.rds")
args <- commandArgs(trailingOnly = TRUE)

baseline.dt <- readRDS(args[1])
intervention.dt <- readRDS(args[2])

bkeys <- key(baseline.dt)
ikeys <- key(intervention.dt)

## perform effectiveness calcs
eff.dt <- intervention.dt[
  # join baseline to interventions on relevant keys
  baseline.dt, on=bkeys, nomatch=0
  # nomatch -> require that results be available for key
  # (e.g. particle & replicate) in both base and intervention
  # i.e., equi-join
][,
  # in data.table joins dt[i.dt, ...], cols in i.dt w/ same name (and not join keys)
  # are prefixed with i.
  .(
    bcases = i.s, icases = s,
    c.bcases = i.c.s, c.icases = c.s,
    # handles i.s == s == 0 case => 1.0 eff
    # cases with i.s == 0, s != 0 assigned Inf
    eff = ifelse(i.s == s, 1.0, (i.s-s)/i.s),
    c.eff = ifelse(i.c.s == c.s, 1.0, (i.c.s-c.s)/i.c.s),
    # handles i.s == s == 0 case => 1 multiplier => 0 log-multiplier
    # cases with i.s == 0, s != 0 assigned Inf
    logm = log(ifelse(i.s == s, 1.0, s/i.s)),
    c.logm = log(ifelse(i.c.s == c.s, 1.0, c.s/i.c.s))
  ),
  keyby=ikeys
]

saveRDS(eff.dt, tail(args, 1))