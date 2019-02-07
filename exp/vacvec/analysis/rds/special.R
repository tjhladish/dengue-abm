# parses the sqlite databases that result from HPC runs
# uses definitions set in projref.R (proJECT refERENCE)

# make command line output less verbose
suppressPackageStartupMessages({
  require(data.table)
  require(RSQLite)
})

# developer args
args <- c("../utils.R", '~/Dropbox/who/mpeak_incidence_response-tmp.sqlite', "~/Dropbox/who/mpeak_intros-tmp.out")

# actual args when used with shell
args <- commandArgs(trailingOnly = TRUE)

# TODO just read in SQL files

# db should be 2nd-to-last arg

qry <- "SELECT M.*, foi
  FROM met M JOIN par P USING(serial) JOIN job J USING(serial)
  WHERE status == 'D';"

source(args[1])

dt <- dbutil(args[2], qry)
dropcols <- grep("_",names(dt))
if (length(dropcols)) dt <- dt[,.SD,.SDcols=-dropcols]
mlt <- melt.data.table(dt, id.vars = c("serial", "foi"))

mlt[, measure := gsub("\\d+","",variable) ]
mlt[, year := as.integer(gsub("^[si]","",variable)) ]
mlt$variable <- NULL


rdt <- dcast.data.table(mlt, serial + year + foi ~ measure, value.var = "value")

intros <- fread(args[3])[, year := year - 101L][year >= 0]
plot.dt <- rdt[
  intros,
  on=.(serial, year), nomatch = 0
]

setnames(plot.dt, "introduced_infections", "intro.i")

plot.dt[, intro.s := floor((s/i)*intro.i) ]
plot.dt[, local.i := i - intro.i ]
plot.dt[, local.s := s - intro.s ]

res <- melt.data.table(plot.dt, id.vars = c("serial","year","foi"))[,
  dtquantiles(probs = c(.25,.5,.75), value, c("lo","med", "hi")),
  by=.(foi, variable, year)
]

tar <- tail(args, 1)

saveRDS(res, tar)