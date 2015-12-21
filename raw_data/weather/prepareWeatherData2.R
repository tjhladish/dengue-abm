# prepare weather data sets
# want two data sets, of daily (TMAX, TMIN) on Date (Date obj, but YYYY-MM-DD format)
# also want day of year, w/ leap days discarded

rm(list=ls())

require(ggplot2)
require(reshape2)
require(data.table)
require(lubridate)

readIn <- function(pth) setkey(fread(
  pth,
  select = c("Year Local", "Month Local", "Day Local", "Hour Local","Temperature (C)"),
  col.names = c('year','month','day','hour','celsius')
)[,
  date := dmy(paste(day,month,year))+hm(paste(hour,"00"))
][, 
  doy := yday(as.Date(paste("2001",month,day,sep="/"))) - ifelse(month > 2, 1, 0)
][!(day == 29 & month == 2)], year, month, day, hour, date, doy)

merida <- readIn("~/Downloads/Hourly/AllYears/32569_Manuel_Crescencio_Rejón_International_Airport__Hourly_1948_2015.csv")
miami <- readIn("~/Downloads/Hourly/AllYears/30883_Miami_International_Airport__Hourly_1948_2015.csv")

firstday <- merida[!is.na(celsius), min(date)]
lastday <- merida[!is.na(celsius), max(date)]

joint <- merida[between(date, firstday, lastday)][miami[between(date, firstday, lastday)]]
t.model <- lm(celsius ~ i.celsius, data=joint)
censor <- as.integer(names(which(abs(rstandard(t.model)) > 5)))
joint[,state:=factor('original',levels=c('original','missing','outlier'))]
joint[is.na(celsius), state := factor('missing',levels=c('original','missing','outlier'))]
orig_outliers <- joint[censor]
joint[censor, celsius := NA][censor, state:= factor('outlier',levels=c('original','missing','outlier'))]
t.model2 <- lm(celsius ~ i.celsius, data=joint)
miami.missing <- joint[is.na(celsius) & is.na(i.celsius), which=TRUE]
miami.interpolate <- (
  joint[miami.missing-1, i.celsius] +
  joint[miami.missing+1, i.celsius]
)/2
joint[miami.missing, i.celsius := miami.interpolate]
joint[is.na(celsius), celsius := predict(t.model2, newdata = .SD)]

dateRange <- list(start=as.POSIXct("1979/1/1 00:00:00", tz = "UTC")-1, end=as.POSIXct("2014/1/1 00:00:00",tz="UTC")-1)

## at this point, slight differences in reconstructed values, due to initial (here) vs later (in orig) pruning by date

# from Chan-Johansson
beta0 = 2.9
betaT = -0.08
vr = 1/4.9

eipfun = function(tempC, beta0, betaT, vr) exp(exp(beta0 + betaT*tempC) + vr/2)

joint[,
  EIP := eipfun(celsius, beta0, betaT, vr)
][,
  EIR := 1/EIP
][,
  EIprogress := EIR/24
][,
  globalindex := .I
]

limdate = joint[dim(joint)[1]:1,list(datelim=cumsum(EIprogress), date)][datelim > 1][1,date]-1

accum <- function(here, expectedMax=30) with(
  joint[
    globalindex >= here
  ][
    1:(24*expectedMax),
    list(cumeir=cumsum(EIprogress), index=.I, add=EIprogress)
  ][
    cumeir > 1,
    list(index, cumeir, add)
  ][1],
  (index - (1-(cumeir - 1)/add))/24
)

lookaheads = joint[
  between(date,dateRange$start, min(limdate,dateRange$end))
][,
  accum(globalindex, 20),
  keyby=list(year, doy, hour)
]

ggplot(lookaheads[, mean(V1), keyby=doy]) + aes(x=doy, y=V1) + geom_line()

# todo this better:
#  - start at several positions (based on number of available threads)
#  - for a position, calculate its required lookahead
#  - next position, has at least that lookahead-1 => can reduce calculations to starting there

# ggplot(joint) + aes(x=date, y=celsius, color=state) + geom_point()