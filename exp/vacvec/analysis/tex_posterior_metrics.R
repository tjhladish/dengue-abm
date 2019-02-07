suppressPackageStartupMessages({
  require(data.table)
  require(jsonlite)
})

args <- c("utils.R", "rds/posterior.rds", "metkey.json", "tex/posterior_met.tex")
args <- commandArgs(trailingOnly = TRUE)

source(args[1])

mqs <- readRDS(args[2])$mets
renamekey <- read_json(args[3])

mqs[, metric := renamekey[[as.character(variable)]], by=variable ]

rws <- mqs[, .(rows=sprintf("%s & %3g & %3g & %3g & [%3g,%3g]",
  metric,
  signif(observed, 3),
  signif(mean, 3),
  signif(med, 3),
  signif(lo, 3),
  signif(hi, 3)
)), by=variable ]

out <- tabularxinnards(rws)

cat(
  out,
  file=tail(args, 1)
)
