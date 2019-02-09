# parses the sqlite databases that result from HPC runs
# uses definitions set in projref.R (proJECT refERENCE)

# make command line output less verbose
suppressPackageStartupMessages({
  require(data.table)
  require(RSQLite)
})

# developer args
args <- c("../utils.R", "projref.rda", "~/Dropbox/who/vacvec-new_yuc-ivn_lag258-reduced_catchup.sqlite", "lag_intervention.rds")
args <- c("../utils.R", "projref.rda", '~/Dropbox/VGRAND/vacvec-new_yuc.sqlite', "baseline.rds")
args <- c("../utils.R", "projref.rda", '~/Dropbox/VGRAND/vacvec-new_yuc-alt_foi.sqlite', "foi_intervention.rds")

# actual args when used with shell
args <- commandArgs(trailingOnly = TRUE)

# TODO turn into generic cabp-pkg
source(args[1])

# bring in reference project definitions
load(args[2])

# TODO just read in SQL files

# db should be 2nd-to-last arg
srcdb <- tail(args, 2)[1]

# target rds should be last arg
# will sniff target name to determine proper keys
tar <- tail(args, 1)





################################# GET THE DATA #################################

# TODO change to input dependency
# build up SQLite query based on which results being parsed
# samplecols from projref.R: c(particle, replicate)
selcols <- c(
  "M.*", # TODO could be more parsimonious here to speed up, e.g., digesting interventions
  sprintf(c("posterior AS %s", "CAST(realization AS INT) AS %s"), samplecols)
)
filters <- c("status == 'D'")

## applies to basic, foi, lag (intervention only)
if (grepl("baseline", tar)) {
  # need no additional columns
  # want only results with no vector control AND no vaccine
  filters <- c(filters, "vector_control == 0", "vac == 0")
} else if (grepl("intervention", tar)) {
  # need the scenario columns
  selcols <- c(selcols, "vector_control AS vc", "vac", "vc_coverage*100 AS vc_coverage", "vac_mech", "catchup")
  # want only results with some intervention
  filters <- c(filters, "(vector_control == 1 OR vac == 1)")
}

## only special analyses foi, lag
if (grepl("foi", tar)) {
  # need extra scenario column
  # applies to both baseline & intervention
  selcols <- c(selcols, "foi")
} else if (grepl("lag", tar)) {
  # need extra scenario column
  selcols <- c(selcols, "vac_first", "ivn_lag")
}

## assemble pieces into query
qry <- sprintf(
  "SELECT %s FROM %s WHERE %s;", # select COLS from TABLE+JOINS where FILTER
  paste(selcols, collapse=", "),
  "met M JOIN par P USING(serial) JOIN job J USING(serial)",
  paste(filters, collapse=" AND ")
)

# dbutil from utils.R
tar.dt <- dbutil(srcdb, qry)





################################# NORMALIZATION #################################

# columns associated with pre-intervention data & job serial
rmv <- c(grep("s_|imm\\d__", names(tar.dt), value = T), "serial")
# data.table syntax for drop columns
tar.dt <- tar.dt[,.SD,.SDcols=-rmv]

# from projref.R
keycols <- samplecols

# if parsing intervention file, add keys + translate scenario info
# trans_% functions from projref.R
if (grepl("intervention", tar)) {
  keycols <- c("scenario", "vc_coverage","vaccine","catchup", keycols)
  tar.dt[vc == 0, vc_coverage := 0]
  # translate vaccine & catchup to standardized lingo for post processing
  trans_vaccine.data.table(tar.dt)
  trans_catchup.data.table(tar.dt)
  trans_scnario.data.table(tar.dt)
  tar.dt$vac_mech <- NULL
  tar.dt$vac <- NULL
  tar.dt$vc <- NULL
}

# add extra keys for special analyses
if (grepl("foi", tar)) keycols <- c("foi", keycols)
if (grepl("lag", tar)) keycols <- c("vac_first", "ivn_lag", keycols)






############################## TO LONG(ISH) FORMAT #############################

# when the data.table gets melted to long format
# will have a variable column, with values corresponding to wide cols
# (s_00, s_01, ...; imm0_00, imm0_01, ...) names
# need to extract the measure (before _) and year (after _)
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
tar.mlt <- parse.meas.yr(
  melt.data.table(tar.dt, id.vars = keycols, variable.factor = FALSE)
)

# recast the data.table to (keys, year, imm proportions, s) cols
result.dt <- dcast.data.table(tar.mlt,
  as.formula(paste(
    paste(c(keycols, "year"), collapse=" + "),
    "measure", sep = " ~ ")
  ),
  value.var = "value"
)

# rename imm0 to seronegative
setnames(result.dt, "imm0", "seronegative")




############################ FINAL CALC & SAVE #########################################

# compute cumulative incidence for each scenario
result.dt[order(year),
  c.s := cumsum(s),
  by = keycols
]

keycols <- c(keycols, "year")
setkeyv(result.dt, keycols)

saveRDS(result.dt, tar)
