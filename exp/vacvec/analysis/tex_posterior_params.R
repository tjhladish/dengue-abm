suppressPackageStartupMessages({
  require(data.table)
  require(knitr)
})

args <- c("rds/posterior.rds", "tex/posterior_par.tex")
args <- commandArgs(trailingOnly = TRUE)

pqs <- readRDS(args[1])$par

cat(
  kable(pqs[,.(`Median (25%,75%)`=sprintf("%.3g (%.3g,%.3g)",med,lo,hi)),
            by=.(Parameter=variable)], format="latex"),
  file=tail(args, 1)
)