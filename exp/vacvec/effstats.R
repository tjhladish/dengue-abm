suppressPackageStartupMessages({
  require(data.table)
  source("utils.R")
  source("projref.R")
})

# debugging args for interactive use
args <- c("comboeff.rds", "effstats.rds")
args <- c("foi_comboeff.rds", "foi_effstats.rds")
args <- c("lag_comboeff.rds", "lag_effstats.rds")
args <- commandArgs(trailingOnly = TRUE)

# read in the relevant comboeff.rds file
comboeff.dt <- readRDS(args[1])
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

saveRDS(res, args[2])