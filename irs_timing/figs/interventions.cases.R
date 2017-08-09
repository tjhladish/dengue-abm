
args <- commandArgs(trailingOnly = TRUE)
# args stopping-interventions.rds interventions.cases.rds

require(data.table)

interventions.dt <- setkey(
  readRDS(args[1])[end_year==50], ## not looking at stopping interventions here
  coverage, particle, year
)

res <- interventions.dt[,
  list(cases=s, cum.cases=cumsum(s), year), by=.(coverage, particle)
]

saveRDS(res, args[2])