suppressPackageStartupMessages({
  require("data.table")
})

args <- commandArgs(trailingOnly = TRUE)

ref <- readRDS(args[1])[vaccine == "edv" & catchup == "vc+vac" & vc_coverage == 75 & variable == "combo.eff"]
cat("ref[year < 10, min(med)]:\n")
cat(ref[year < 10, min(med)],"\n")
cat("ref[year < 20, min(med)]:\n")
cat(ref[year < 20, min(med)],"\n")
cat("ref[year < 40, min(med)]:\n")
cat(ref[year < 40, min(med)],"\n")