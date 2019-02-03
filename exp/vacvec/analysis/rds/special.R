# parses the sqlite databases that result from HPC runs
# uses definitions set in projref.R (proJECT refERENCE)

# make command line output less verbose
suppressPackageStartupMessages({
  require(data.table)
  require(RSQLite)
})

# developer args
args <- c("../utils.R", '~/Dropbox/who/mpeak_incidence_response.sqlite', "~/Dropbox/who/mpeak_intros.out")

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
dt <- dt[,.SD,.SDcols=-dropcols]
mlt <- melt.data.table(dt, id.vars = c("serial", "foi"))

mlt[, measure := gsub("\\d+","",variable) ]
mlt[, year := as.integer(gsub("^[si]","",variable)) ]
mlt$variable <- NULL
rdt <- dcast.data.table(mlt, serial + year + foi ~ measure, value.var = "value")

intros <- fread(args[3])
plot.dt <- rdt[
  intros[, year := year - 101L][year >= 0],
  on=.(serial, year), nomatch = 0
]

it.dt <- plot.dt[, .(y=as.numeric(median(i, na.rm = TRUE)), measure="infections", fraction="total"), keyby=.(foi, year)]
st.dt <- plot.dt[, .(y=as.numeric(median(s, na.rm = TRUE)), measure="cases", fraction="total"), keyby=.(foi, year)]
m.dt <- plot.dt[, .(m=median(s/i, na.rm = TRUE)), keyby=.(foi,year)]

ii.dt <- plot.dt[, .(y=as.numeric(median(introduced_infections, na.rm = TRUE)), measure="infections", fraction="introduced"), keyby=.(foi, year)]
si.dt <- ii.dt[m.dt, on=.(foi,year)][,.(y=y*m, measure="cases", fraction), keyby=.(foi, year)]

bind.dt <- rbind(it.dt, st.dt, ii.dt, si.dt)

tar <- tail(args, 1)

saveRDS(bind.dt, tar)