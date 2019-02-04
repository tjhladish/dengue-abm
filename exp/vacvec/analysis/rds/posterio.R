# parses the sqlite databases that result from HPC runs
# uses definitions set in projref.R (proJECT refERENCE)

# make command line output less verbose
suppressPackageStartupMessages({
  require(data.table)
  require(RSQLite)
})

# developer args
args <- c("../utils.R", '~/Dropbox/who/new_yuc_posterior.sqlite', "posterior.rds")

# actual args when used with shell
args <- commandArgs(trailingOnly = TRUE)

# TODO just read in SQL files

# db should be 2nd-to-last arg

mqry <- "SELECT M.*
  FROM met M JOIN job J USING(serial)
  WHERE status == 'D';"

pqry <- "SELECT U.*
  FROM upar U JOIN job J USING(serial)
  WHERE status == 'D';"


source(args[1])

mdt <- melt.data.table(dbutil(args[2], mqry), id.vars = "serial")
pdt <- melt.data.table(dbutil(args[2], pqry)[, seed := NULL ], id.vars = "serial")

mqs <- mdt[,dtquantiles(c(.25,.5,.75), value, pnames = c("lo","med","hi")), by=variable]
pqs <- pdt[,dtquantiles(c(.25,.5,.75), value, pnames = c("lo","med","hi")), by=variable]

saveRDS(list(pars=pqs,mets=mqs), tail(args, 1))