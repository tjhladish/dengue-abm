suppressPackageStartupMessages({
  require(data.table)
  require(jsonlite)
})

args <- c("utils.R", "rds/posterior.rds", "parkey.json", "tex/posterior_par.tex")
args <- commandArgs(trailingOnly = TRUE)

source(args[1])

pqs <- readRDS(args[2])$pars
renamekey <- read_json(args[3])

pqs[, parameter := renamekey[[as.character(variable)]], by=variable ]

rws <- pqs[, .(rows=sprintf("%s & %s & %s & %.3g & [%.3g,%.3g]",#    %3g & %3g & %3g & $[%3g,%3g]$
  parameter,
  ifelse(distro=="NORMAL",
    sprintf("$\\mathcal{N}$(%.7g, %.7g)",distro_par1,distro_par2),
    "?"
  ),
  ifelse(uxf != "IDENTITY",
    sprintf("%.3g & [%.3g, %.3g]", uxfmean, uxflower, uxfupper),
    " & "
  ),
  signif(med, 3),
  signif(lo, 3),
  signif(hi, 3)
)), by=variable ]

out <- tabularxinnards(rws)

cat(
  out,
  file=tail(args, 1)
)
