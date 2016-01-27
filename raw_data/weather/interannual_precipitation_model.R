all_data = read.csv("NOAA_yucatan_daily.csv")

# Filter on data from the Merida airport
d = all_data[all_data$NAME=='AEROP.INTERNACIONAL',]
d$DATE = as.Date(d$DATE)

# lag by one week to account for egg-to-emergence time, based on Table 2 in Rueda et al, Temperature-dependent
# development and survival rates of Culex quinquefasciatus and Aedes aegypti (Diptera: Culicidae), J Med Ent, 1990, 27:5
mosquito_lag = 7
# Filter on dates of interest
precip=data.frame(DATE=seq(as.Date("1978/4/10")-mosquito_lag, as.Date("2014/4/9")-mosquito_lag, "days"))
precip = precip[-which(strftime(precip$DATE, '%m-%d') == '02-29'),,drop=FALSE]
precip$PRCP = d$PRCP[match(precip$DATE, d$DATE)]

precip$bool0 = precip$PRCP > 0
precip$month = strftime(precip$DATE, '%m')
precip$day = strftime(precip$DATE, '%m-%d')
precip$year = strftime(precip$DATE, '%Y')
num_epiyears = length(1978:2013) # starts apr 10, 1978, ends apr 9, 2014
precip$epiyear = rep(1:num_epiyears, each=365)

# aggregate by day of year (fraction of years with rain on each day)
seasonal_rain = aggregate(precip$bool0, by=list(precip$day), mean, na.rm=T)

# aggregate by year (fraction of days with rain in each year)
yearly_multipliers = aggregate(precip$bool0, by=list(precip$epiyear), mean, na.rm=T)$x
early_mean = mean(yearly_multipliers[1:10])
late_mean = mean(rev(yearly_multipliers)[1:10])
# pad sequence so smoothing function predicts reasonably past the end points of the known series
yearly_multipliers = rep(c(early_mean, yearly_multipliers, late_mean), each=365)

names(seasonal_rain) = c('day', 'precip')

num_years = 3 # pad an extra year at beginning and end, so that we can truncate and wrap
seasonal_tmp = data.frame(precip=rep(early_mean*seasonal_rain$precip, num_years), index=1:(num_years*365))
loess_model.burnin_period = loess(precip~index, data=seasonal_tmp, span=0.3/num_years, degree=2, evaluation=round(50*num_years))
smooth.burnin = predict(loess_model.burnin_period, newdata = 1:365 + 365 - mosquito_lag)
normal_max = max(smooth.burnin)
write.table(data.frame(smooth=smooth.burnin/normal_max, day=1:365),
            'burnin_precipitation.txt', row.names=F, col.names=F, quote=F)


# Epidemic year is assumed to start/end on day 100/99, so we need to offest the precip data
# This has the added benefit of putting transitions at the low point in the precip curve, thus minimizing
# transition artifacts
offset_seasonal_rain = rbind(seasonal_rain[-c(1:99),], seasonal_rain[c(1:99),])
# plus 2 on following line because we pad the beginning and end with a year that is later discarded after smoothing
sr_copy = data.frame(day=rep(offset_seasonal_rain$day, num_epiyears+2), precip=rep(offset_seasonal_rain$precip, num_epiyears+2))
sr_copy$daily = sr_copy$precip*yearly_multipliers
sr_copy$index = 1:length(sr_copy$daily)
nrow = dim(sr_copy)[1]

num_years = length(sr_copy$daily)/365
loess_model.fitting_period = loess(daily~index, data=sr_copy, span=0.3/num_years, degree=2, evaluation=round(50*num_years))
smooth.entire = predict(loess_model.fitting_period, newdata = sr_copy$index)
burnin = 1:(265 + 365)
good_days = 1:(365*35+100)
smooth.cropped = smooth.entire[-burnin]
sr_copy = sr_copy[-burnin,]
smooth.cropped = smooth.cropped[good_days]
sr_copy = sr_copy[good_days,]

year_vals = rep(1979:2014, each=365)[1:(365*(num_epiyears-1)+100)]
day_vals  = rep(1:365, num_epiyears)[1:(365*(num_epiyears-1)+100)]
print(normal_max)
write.table(data.frame(smooth=smooth.cropped/normal_max, year=year_vals, day=day_vals),
#write.table(data.frame(smooth=smooth.cropped, year=year_vals, day=day_vals),
            'series_precipitation.txt', row.names=F, col.names=F, quote=F)

pdf('rain_refit2.pdf', width=20, height=6)
plot(sr_copy$daily, type='l', lwd=0.5)#, xlim=c(4000,5000))
lines(smooth.cropped, col=2, lwd=2)
#lines(seasonal_rain$smooth,col=2, lwd=2)
#lines(1:365+8, seasonal_rain$smooth,col=3, lwd=2)
dev.off()
