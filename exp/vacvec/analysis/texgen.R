suppressPackageStartupMessages({
  require(data.table)
})

.args <- c("defs.template", "rds/effstats.rds", "rds/lag_effstats.rds", "numerical_results.tex")
.args <- commandArgs(trailingOnly = TRUE)

template <- paste(readLines(.args[1]), collapse="\n")

ref.dt <- readRDS(.args[2])
lagref.dt <- readRDS(.args[3])

ideal.slice <- ref.dt[
  year == 9 & variable == "combo.eff" &
  vaccine == "edv" & vc_coverage == 75 &
  catchup == "vc+vac"
]

lag.slice <- lagref.dt[variable == "combo.eff"][]

res <- sprintf(template,
  100*ideal.slice$med, 100*ideal.slice$lo, 100*ideal.slice$hi
)

cat(res, file = tail(.args,1))