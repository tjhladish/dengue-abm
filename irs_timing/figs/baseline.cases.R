
args <- commandArgs(trailingOnly = TRUE)

require(data.table)

baseline.dt <- setkey(readRDS(args[1]), particle, year)

res <- baseline.dt[,
  list(cases=s, cum.cases=cumsum(s), year), by=particle
]

saveRDS(res, args[2])