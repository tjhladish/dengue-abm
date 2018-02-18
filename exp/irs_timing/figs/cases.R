
args <- commandArgs(trailingOnly = TRUE)
# args <- c("~/Dropbox/who/fig1_data/stopping-interventions.rds", "~/Dropbox/who/fig1_data/interventions.cases.rds")
# args <- c("~/Dropbox/who/fig1_data/stopping-baseline.rds", "~/Dropbox/who/fig1_data/baseline.cases.rds")
require(data.table)

ref.dt <- readRDS(args[1])
# assert: has (year, s) columns (corresponding to x, y cols), + other key cols
keys <- grep("^(imm|s)",names(ref.dt), invert = T, value = T)
setkeyv(ref.dt, keys)
notyear <- grep("^year$", keys, invert = T, value = T)
# need to ensure rows are sorted by year

# interventions.dt <- setkey(
#   readRDS(args[1])[end_year==50], ## not looking at stopping interventions here
#   coverage, particle, year
# )

res <- ref.dt[,
  list(cases=s, cum.cases=cumsum(s), year), by=notyear
]

saveRDS(res, args[2])