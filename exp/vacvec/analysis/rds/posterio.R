# parses the sqlite databases that result from HPC runs
# uses definitions set in projref.R (proJECT refERENCE)

# make command line output less verbose
suppressPackageStartupMessages({
  require(data.table)
  require(RSQLite)
  require(jsonlite)
})

# developer args
args <- c("../utils.R", '~/Dropbox/who/new_yuc_posterior.sqlite', "../../../abc-irs_refit2/abc-irs_refit2.json", "posterior.rds")

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

mqs <- mdt[,dtquantiles(c(.025,.5,.975), value, pnames = c("lo","med","hi"), mn=T), by=variable]
pqs <- pdt[,dtquantiles(c(.025,.5,.975), value, pnames = c("lo","med","hi"), mn=T), by=variable]

refjson <- read_json(args[3])

metref <- refjson$metrics
vals <- lapply(metref, function(i) i$value)
names(vals) <- sapply(metref, function(i) i$name)
mqs[, observed := as.numeric(vals[[as.character(variable)]]), by=variable]

parref <- refjson$parameters
parnames <- sapply(parref, function(i) ifelse(is.null(i$short_name),i$name,i$short_name))
names(parref) <- parnames

pqs[, distro := parref[[as.character(variable)]]$dist_type, by=variable ]
pqs[, distro_par1 := as.numeric(parref[[as.character(variable)]]$par1), by=variable ]
pqs[, distro_par2 := as.numeric(parref[[as.character(variable)]]$par2), by=variable ]
pqs[, uxf := ifelse(
  is.null(parref[[as.character(variable)]]$untransform),
  "IDENTITY",
  ifelse(is.character(parref[[as.character(variable)]]$untransform),
         parref[[as.character(variable)]]$untransform,
         parref[[as.character(variable)]]$untransform$type
  )
), by=variable ]

pqs[, uxflower := ifelse(
  is.null(parref[[as.character(variable)]]$untransform),
  NA_real_,
  as.numeric(ifelse(is.character(parref[[as.character(variable)]]$untransform),
         0,
         parref[[as.character(variable)]]$untransform$min
  ))
), by=variable ]
pqs[, uxfupper := ifelse(
  is.null(parref[[as.character(variable)]]$untransform),
  NA_real_,
  as.numeric(ifelse(is.character(parref[[as.character(variable)]]$untransform),
         1,
         parref[[as.character(variable)]]$untransform$max
  ))
), by=variable ]

pqs[uxf == "LOGISTIC" & distro == "NORMAL", uxfmean := logistic(distro_par1, L=uxfupper-uxflower) + uxflower]

saveRDS(list(pars=pqs,mets=mqs), tail(args, 1))