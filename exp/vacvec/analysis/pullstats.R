suppressPackageStartupMessages({
  require("data.table")
})

args <- c("rds/effstats.rds", "rds/lag_effstats.rds")
args <- commandArgs(trailingOnly = TRUE)

ref <- readRDS(args[1])

slice1 <- ref[vaccine == "edv" & catchup == "vc+vac" & vc_coverage == 75 & variable == "combo.eff"]

cat("edv w/ catchup, 75% vc, years 0:9, min(median eff):\n")
cat(slice1[year < 10, min(med)],"\n")
cat("edv w/ catchup, 75% vc, years 0:19, min(median eff):\n")
cat(slice1[year < 20, min(med)],"\n")
cat("edv w/ catchup, 75% vc, years 0:39, min(median eff):\n")
cat(slice1[year < 40, min(med)],"\n")

slice2 <- ref[vaccine == "edv" & catchup == "vc+vac" & vc_coverage == 75 & variable == "c.combo.eff" & year == 9]

cat("edv w/ catchup, 75% vc, year 9 cumulative eff (95% PI):\n")
cat(slice2[, c(med, lo.lo, hi.hi)],"\n")

ref2 <- readRDS(args[2])[variable == "c.combo.eff" & year == 19]

cat("lag interventions, min c.combo.eff at 19 years:\n")
cat(ref2[, min(med)])