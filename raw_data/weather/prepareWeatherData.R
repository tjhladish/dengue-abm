# prepare weather data sets
# want two data sets, of daily (TMAX, TMIN) on Date (Date obj, but YYYY-MM-DD format)
# also want day of year, w/ leap days discarded

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

# check weather data visually

require(ggplot2)
require(reshape2)

both <- melt(rbind(
  cbind(merida, location=factor("MERIDA", levels=c("MERIDA","MIAMI"), ordered = T)),
  cbind(miami, location=factor("MIAMI", levels=c("MERIDA","MIAMI"), ordered = T))
), id.vars = c("doy","location","DATE"), value.name = "celsius", variable.name = "extrema", na.rm = T)

ggplot(both) + theme_bw() + facet_grid(. ~ location, scales = "free") + 
  aes(x=doy, y=celsius, color=extrema, group=interaction(year(DATE), extrema)) + geom_line(alpha=0.1) +
  scale_color_manual(values=c(TMAX='red',TMIN='blue'))

