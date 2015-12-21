# prepare weather data sets
# want two data sets, of daily (TMAX, TMIN) on Date (Date obj, but YYYY-MM-DD format)
# also want day of year, w/ leap days discarded

require(ggplot2)
require(reshape2)
require(data.table)


miami <- setkey(fread("NOAA_miami.csv")[!(DAY == 29 & MONTH == 2),
  list(
    TMAX = TMAX/10.0,
    TMIN = TMIN/10.0,
    DATE = as.Date(paste(YEAR,MONTH,DAY,sep="/")),
    doy = yday(as.Date(paste("2001",MONTH,DAY,sep="/")))
  )
], DATE)[!(is.na(TMAX) & is.na(TMIN))]

saveRDS(miami, "miami.RData")

merida = setkey(fread("NOAA_yucatan_daily.csv")[NAME=='AEROP.INTERNACIONAL'][grep("02/29", DATE, invert = T),
  list(
    DATE = as.Date(DATE),
    TMAX = TMAX/10,
    TMIN = TMIN/10,
    doy = yday(as.Date(paste("2001",month(DATE),mday(DATE),sep="/")))
  )
], DATE)[!(is.na(TMAX) & is.na(TMIN))]

saveRDS(merida, "merida.RData")
