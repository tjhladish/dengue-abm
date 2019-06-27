suppressPackageStartupMessages({
  require(data.table)
  require(RSQLite)
})
# one script to merge them all

# FOISIDEDB ?= ~/Dropbox/who/dengue/mpeak_incidence_response-expanded.sqlite
# FOISIDETGZ ?= ~/Dropbox/who/dengue/mpeak_intros-expanded.out
# ABCDB ?= ~/Dropbox/who/dengue/new_yuc_posterior.sqlite

.args <- c(paste0("~/Dropbox/who/dengue/",c(
  "vacvec-new_yuc.sqlite", # original results - no cyd-tdv testing, foi = 1, lag = 0
  "vacvec-new_yuc-alt_foi.sqlite", # foi results - no cyd-tdv testing, lag = 0
  "vacvec-new_yuc-ivn_lag258-reduced_catchup.sqlite", # lag results; no cyd-tdv (so testing = NA); foi = 1
  "vacvec-new_yuc-CYD_TDV-seropos_only.sqlite", # cyd-tdv with perfect testing; foi = 1
  "vacvec-new_yuc-alt_foi-CYD_TDV-seropos_only.sqlite" # cyd-tdv with perfect testing, different fois
)), "~/Dropbox/cabpshare/vacvec-CYD_TDV-w_test.sqlite")

.args <- c(.args, "baseline.rds") # select only non-intervention items
.args <- c(.args, "intervention.rds") # select only intervention items

.args <- commandArgs(trailingOnly = T)

tar <- tail(.args, 1)

samplecols <- c("particle", "replicate")
dbutil <- function(dbfile, sql) {
  drv = dbDriver("SQLite")
  db = dbConnect(drv, dbfile, flags=SQLITE_RO)
  res <- data.table(dbGetQuery(db, sql))
  dbDisconnect(db)
  return(res)
}

# sum(unlist(lapply(.args[1:6], function(dbf) dbGetQuery(dbConnect(dbDriver("SQLite"), dbf),"SELECT COUNT(*) FROM met JOIN job USING(serial) WHERE status=='D';"))))

selcols <- c(
  paste(sprintf("s%02i",0:39), collapse = ", "),
  sprintf(c("posterior AS %s", "CAST(realization AS INT) AS %s"), samplecols)
)
filters <- c("status == 'D'")

is_baseline <- grepl("baseline", tar)

keycols <- c(samplecols, "foi")

if (is_baseline) {
  # need no additional columns
  # want only results with no vector control AND no vaccine
  filters <- c(filters, "vector_control == 0", "vac == 0")
} else {
  # need the scenario columns
  selcols <- c(selcols, "vector_control AS vc", "vac", "vc_coverage*100 AS vc_coverage", "vac_mech", "catchup")
  # want only results with some intervention
  filters <- c(filters, "(vector_control == 1 OR vac == 1)")
}

substs <- c(
  nolag="0 AS ivn_lag, NULL AS vac_first",
  lag="ivn_lag, vac_first",
  nofoi="1.0 AS foi",
  foi="foi",
  notest="NULL AS false_pos, NULL AS false_neg",
  ptest="0 AS false_pos, 0 AS false_neg",
  rtest="false_pos, false_neg"
)

qry <- function(scols, fs) sprintf(
  "SELECT %s FROM %s WHERE %s;", # select COLS from TABLE+JOINS where FILTER
  paste(scols, collapse=", "),
  "met M JOIN par P USING(serial) JOIN job J USING(serial)",
  paste(fs, collapse=" AND ")
)

queries <- character(6)

notest_nofoi_nolag <- .args[1]
queries[1] <- { if (is_baseline) { 
  qry(c(selcols, substs["nofoi"]),filters)
} else {
  qry(c(selcols, substs[c("nofoi", "nolag", "notest")]),filters)
}}

notest_foi_nolag <- .args[2]
queries[2] <- { if (is_baseline) { 
  qry(c(selcols, substs["foi"]), filters)
} else {
  qry(c(selcols, substs[c("foi", "nolag", "notest")]),filters)
}}

notest_nofoi_lag <- .args[3]
queries[3] <- { if (is_baseline) { 
  qry(c(selcols, substs["nofoi"]), filters)
} else {
  qry(c(selcols, substs[c("nofoi", "lag", "notest")]),filters)
}}

ptest_nofoi_nolag <- .args[4]
queries[4] <- { if (is_baseline) { 
  qry(c(selcols, substs["nofoi"]), filters)
} else {
  qry(c(selcols, substs[c("nofoi", "nolag", "ptest")]),filters)
}}

ptest_foi_nolag <- .args[5]
queries[5] <- { if (is_baseline) { 
  qry(c(selcols, substs["foi"]), filters)
} else {
  qry(c(selcols, substs[c("foi", "nolag", "ptest")]),filters)
}}

rtest_foi_nolag <- .args[6]
queries[6] <- { if (is_baseline) { 
  qry(c(selcols, substs["foi"]), filters)
} else {
  qry(c(selcols, substs[c("foi", "nolag", "rtest")]),filters)
}}

bindthemall <- rbindlist(mapply(dbutil, dbfile=.args[1:6], sql=queries, SIMPLIFY = F))

load("projref.rda")

if (!is_baseline) {
  keycols <- c(keycols, "scenario", "vaccine", "vc_coverage","catchup","ivn_lag","vac_first","false_pos","false_neg")
  bindthemall[vc == 0, vc_coverage := 0]
  # translate vaccine & catchup to standardized lingo for post processing
  trans_vaccine.data.table(bindthemall)
  trans_catchup.data.table(bindthemall)
  trans_scnario.data.table(bindthemall)
  bindthemall$vc <- bindthemall$vac <- bindthemall$vac_mech <- NULL
}


parse.meas.yr <- function(dt) {
  dt[,
     measure := gsub("(s|imm\\d).+","\\1", variable)
     ][,
       year    := as.integer(gsub("(s|imm[0-4]_)","", variable))
       ]
  dt$variable <- NULL
  return(dt)
}

# melt the data.table, then parse it to get measure and year
mlt <- parse.meas.yr(
  melt.data.table(bindthemall, id.vars = keycols, variable.factor = FALSE)
)

# recast the data.table to (keys, year, imm proportions, s) cols
result.dt <- dcast.data.table(mlt,
                              as.formula(paste(
                                paste(c(keycols, "year"), collapse=" + "),
                                "measure", sep = " ~ ")
                              ),
                              value.var = "value"
)

result.dt[order(year),
          c.s := cumsum(s),
          by = keycols
          ]

keycols <- c(keycols, "year")
setkeyv(result.dt, keycols)

saveRDS(result.dt, tar)
