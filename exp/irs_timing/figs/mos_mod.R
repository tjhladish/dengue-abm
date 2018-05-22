require(reshape2)
require(data.table)
require(lubridate)

#args <- c("~/Dropbox/who/fig1_data/mos.rds", "~/Dropbox/who/fig1_data/mod-mos.rds")
args <- commandArgs(trailingOnly = TRUE)

#basemospop <- readRDS(args[1])

makexform <- function(
  start=150L, coverage=0.75,
  efficacy=0.8, duration=90,
  durability=90, householdfrac=376400/475362
) {
  initial <- data.table(
    doy = 1L:365L
  )
  
  if (duration == 365) { # year long rollout - start actually irrelevant
    initial[, multiplier := (1-householdfrac) + householdfrac * coverage * (durability / 365)* efficacy] 
  } else {
    # only care about a certain time frame - start => duration + durability.
    # might wrap, so work with alternate zero?
    # currently ignores the possibility of starting next year treatment while current still active
    
    coverage_increase <- c((1:duration)/duration, rep(1, durability-1))
    coverage_decline <- c(rep(0, durability-1), (1:duration)/duration)
    netcoverage = coverage_increase - coverage_decline
    effect <- netcoverage*coverage*efficacy
    
    add <- data.table(
      doy = (start - 1L + (1L:length(effect))) %% 365L,
      multiplier = ((1-householdfrac) + householdfrac * (1-effect))
    )
    
    add[doy==0, doy := 365L]
    
    initial <- merge(initial, add, on="doy", all = TRUE)
    
    initial[is.na(multiplier), multiplier := 1]
  }
  
  return(initial)
}

scenarios <- data.table(expand.grid(
  start = c(148L, 323L),
  coverage = (1:3)/4,
  efficacy = c(0.4, 0.6, 0.8),
  duration = c(1, 90, 365),
  durability = c(30, 90, 150)
))

res <- scenarios[,
  makexform(start, coverage, efficacy, duration, durability),
  by=.(start, coverage, efficacy, duration, durability)
]

saveRDS(res, tail(args,1))
