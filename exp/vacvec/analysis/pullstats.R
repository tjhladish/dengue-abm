suppressPackageStartupMessages({
  require("data.table")
})

args <- c("rds/effstats.rds", "rds/lag_effstats.rds", "rds/foi_v_mpop.rds")
args <- commandArgs(trailingOnly = TRUE)

ref <- readRDS(args[1])

slice1 <- ref[vaccine == "edv" & catchup == "vc+vac" & vc_coverage == 75 & variable == "combo.eff"]

cat("edv w/ catchup, 75% vc, years 0:9, min(median eff):\n")
tenyearminmed <- slice1[year < 10, min(med)]
cat(tenyearminmed,"\n")

cat("edv w/ catchup, 75% vc, years 0:19, min(median eff):\n")
twentyminmed <- slice1[year < 20, min(med)]
cat(twentyminmed,"\n")

cat("edv w/ catchup, 75% vc, years 0:39, min(median eff):\n")
cat(slice1[year < 40, min(med)],"\n")

slice2 <- ref[vaccine == "edv" & catchup == "vc+vac" & vc_coverage == 75 & variable == "c.combo.eff" & year == 9]

cat("edv w/ catchup, 75% vc, year 9 cumulative eff (95% PI):\n")
cat(slice2[, c(med, lo.lo, hi.hi)],"\n")

ref2 <- readRDS(args[2])[variable == "c.combo.eff" & year == 19]

cat("lag interventions, min c.combo.eff at 19 years:\n")
cat(ref2[, min(med)],"\n")

ref3 <- readRDS(args[3])[, .(val=median(value)), by=.(foi, variable)]
mlt.ref3 <- dcast.data.table(ref3, foi ~ variable, value.var = "val")
targetfoi <- mlt.ref3[i < 2*intro.i, max(foi)]
cat("target foi corresponding to infections < 2x intros:\n")
cat(targetfoi,"\n")

refbaseline <- ref3[foi==1,.(ref = val),by=.(foi, variable)]
joinref3 <- ref3[variable %in% c("i","s")][refbaseline, on=.(variable), nomatch=NULL]

cat("compared to 10 year mark:\n")
cat("measure foi val reducedref\n")
joinref3[val < ref*(1-tenyearminmed),cat(as.character(variable),foi,val,"vs",ref*(1-tenyearminmed),"\n"), by=.(variable, foi)]
cat("compared to 20 year mark:\n")
cat("measure foi val reducedref\n")
joinref3[val < ref*(1-twentyminmed),cat(as.character(variable),foi,val,"vs",ref*(1-twentyminmed),"\n"), by=.(variable, foi)]

