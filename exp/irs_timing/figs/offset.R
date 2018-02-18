
require(data.table)
require(lubridate)

## get script args; for debugging, uncomment section that follows
args <- commandArgs(trailingOnly = TRUE)

## read in digested peaks
peaks.dt <- fread(args[1])

fwrite(peaks.dt[,
  peak.mid.date := ymd(paste0("1970 ",max.doy))+floor(duration/2)+floor(durability/2)
][
  duration != 365,
  .(peak.mid.date=format(unique(peak.mid.date),"%b %d")),
  keyby=.(coverage, duration, durability)
], args[2])
