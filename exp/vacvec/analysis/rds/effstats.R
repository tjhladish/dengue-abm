suppressPackageStartupMessages({
  require(data.table)
})

# debugging args for interactive use
args <- c("../utils.R", "projref.rda", "comboeff.rds", "effstats.rds")
args <- c("../utils.R", "projref.rda", "foi_comboeff.rds", "foi_effstats.rds")
args <- c("../utils.R", "projref.rda", "lag_comboeff.rds", "lag_effstats.rds")
args <- commandArgs(trailingOnly = TRUE)

source(args[1])
load(args[2])

# read in the relevant comboeff.rds file
comboeff.dt <- readRDS(args[3])
# ...and introspect keys
ckeys <- key(comboeff.dt)

# melt into long format for easy quantile'ing
mlt <- melt.data.table(
  comboeff.dt,
  id.vars = ckeys
)

# the new keys exclude the samplecols (what will be quantile'd over)
# and add "variable". This is the default name from melt, corresponds
# to pre-melt column names (that aren't id.vars / keys)
mkeys <- c("variable", setdiff(ckeys, samplecols))

res <- mlt[,
  dtquantiles(quantile_probs, value),
  keyby=mkeys
]

saveRDS(res, tail(args, 1))
