suppressPackageStartupMessages({
  require(data.table)
  require(knitr)
})

args <- c("rds/posterior.rds", "tex/posterior_met.tex")
args <- commandArgs(trailingOnly = TRUE)

mqs <- readRDS(args[1])$mets

cat(
  kable(mqs[,.(`Median (25%,75%)`=sprintf("%.3g (%.3g,%.3g)",med,lo,hi)),
  by=.(Metric=variable)], format="latex"),
  file=tail(args, 1)
)
