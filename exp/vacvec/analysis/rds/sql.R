# parses the sqlite databases that result from HPC runs
# uses definitions set in projref.R (proJECT refERENCE)

# developer args
args <- c("projref.rda", "baseline.sql")
args <- c("projref.rda", "intervention.sql")

# actual args when used with shell
args <- commandArgs(trailingOnly = TRUE)

# sources for shared definitions;
load(args[1])

# target sql should be last arg
# use the target name to determine proper keys
tar <- tail(args, 1)


# build up SQLite query based on which results being parsed
# samplecols from projref.rda: c(particle, replicate)
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
  selcols <- c(selcols, "vac_first")
}

## assemble pieces into query
qry <- sprintf(
  "SELECT %s FROM %s WHERE %s;", # select COLS from TABLE+JOINS where FILTER
  paste(selcols, collapse=", "),
  "met M JOIN par P USING(serial) JOIN job J USING(serial)",
  paste(filters, collapse=" AND ")
)

write(qry, file=tar)
