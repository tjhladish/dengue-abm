
require(data.table)
require(lubridate)

## get script args; for debugging, uncomment section that follows
args <- commandArgs(trailingOnly = TRUE)
# args <- c(
#   paste0("~/Dropbox/who/fig1_data/",
#          c("coverage", "duration", "durability"),
#          ".rds"),
#   "~/Dropbox/who/fig1_data/tab1.csv"
# )

## read in assorted input data
coverage.dt <- readRDS(args[1])[layer=="background"]
duration.dt <- readRDS(args[2])[layer=="background"]
durability.dt <- readRDS(args[3])[layer=="background"]

minmax <- function(dt) {
  dt[,.(
    max=max(value), max.doy=format(as_date(0)+doy[which.max(value)]-1,"%b %d"),
    min=min(value), min.doy=format(as_date(0)+doy[which.min(value)]-1,"%b %d")
  ), by=.(coverage, duration, durability)][,
    ratio := round(max/min, 2)
  ][,
    rmax := round(max,2)
  ][,
    rmin := round(min,2)
  ]
}

fwrite(rbind(minmax(coverage.dt),minmax(duration.dt),minmax(durability.dt)), args[4])
