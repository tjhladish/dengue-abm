all_data = read.csv("NOAA_yucatan_daily.csv")

# work with Merida data only
d = all_data[all_data$NAME=='AEROP.INTERNACIONAL',]
d$DATE = as.Date(d$DATE)

precip=data.frame(DATE=seq(as.Date("1978/4/10"), as.Date("2014/4/9"), "days"))
precip = precip[-which(strftime(precip$DATE, '%m-%d') == '02-29'),,drop=FALSE]
precip$PRCP = d$PRCP[match(precip$DATE, d$DATE)]
#precip = precip[-which(precip$DATE==as.Date("2014/1/1")),]

precip$bool0 = precip$PRCP > 0
precip$month = strftime(precip$DATE, '%m')
precip$day = strftime(precip$DATE, '%m-%d')
precip$year = strftime(precip$DATE, '%Y')
precip$epiyear = rep(1:36, each=365)
seasonal_rain = aggregate(precip$bool0, by=list(precip$day), mean, na.rm=T)
yearly_rain = aggregate(precip$bool0, by=list(precip$epiyear), mean, na.rm=T)
yearly_multipliers = yearly_rain$x
mean_val = mean(yearly_multipliers, na.rm=T)
yearly_multipliers = rep(c(mean_val, yearly_multipliers, mean_val), each=365)

names(seasonal_rain) = c('day', 'precip')
offset_seasonal_rain = rbind(seasonal_rain[-c(1:99),], seasonal_rain[c(1:99),])
sr_copy = data.frame(day=rep(offset_seasonal_rain$day, 38), precip=rep(offset_seasonal_rain$precip, 38))
sr_copy$daily = sr_copy$precip*yearly_multipliers
sr_copy$index = 1:length(sr_copy$daily)
nrow = dim(sr_copy)[1]
ver1 = smooth.spline(rep(seasonal_rain$precip, 22)[1:8000], spar=0.1)$y[1:365+365]
ver2 = smooth.spline(rep(seasonal_rain$precip, 3),spar=0.6)$y[1:365+365]
#seasonal_rain$smooth = mod3$y[1:365 + 365]

year_vals = rep(1979:2014, each=365)[1:(365*35+100)]
day_vals  = rep(1:365, 36)[1:(365*35+100)]
#write.table(data.frame(smooth=smooth.cropped, year=year_vals, day=rep(1:365,36)[1:(365*35+100)]),
#            'series_precipitation.txt', row.names=F, col.names=F, quote=F)


num_years = length(sr_copy$daily)/365
loess_model = loess(daily~index, data=sr_copy, span=0.3/num_years, degree=2, evaluation=round(50*num_years))
smooth.entire = predict(loess_model, newdata = sr_copy$index)
#smooth.early = smooth.spline(sr_copy$daily[c(1:8000)], spar=0.1)$y
#smooth.late = smooth.spline(sr_copy$daily[c((nrow-7999):nrow)], spar=0.1)$y
#smooth.merged = c(smooth.early[1:(nrow/2)], smooth.late[1:8000])
#smooth.merged = c(smooth.early[1:(nrow/2)], smooth.late[(8001-(nrow/2)):8000])
#smooth.cropped= smooth.merged
#smooth.cropped= smooth.entire
burnin = 1:(265 + 365)
good_days = 1:(365*35+100)
smooth.cropped = smooth.entire[-burnin]
sr_copy = sr_copy[-burnin,]
smooth.cropped = smooth.cropped[good_days]
sr_copy = sr_copy[good_days,]
#write.table(data.frame(smooth=smooth.cropped, year=rep(1979:2014, each=365)[1:(365*35+100)], day=rep(1:365,36)[1:(365*35+100)]),
#            'series_precipitation.txt', row.names=F, col.names=F, quote=F)
#write.table(seasonal_rain[, c('smooth', 'precip', 'day')], 'merida_preciptiation.out', row.names=F, col.names=F, quote=F)

#seasonal_rain$mosquitoes = seasonal_rain$smooth/max(seasonal_rain$smooth)
# lag by one week to account for egg-to-emergence time, based on Table 2 in Rueda et al, Temperature-dependent
# development and survival rates of Culex quinquefasciatus and Aedes aegypti (Diptera: Culicidae), J Med Ent, 1990, 27:5
#seasonal_rain$mosquitoes = seasonal_rain$mosquitoes[(0:364-7)%%365 + 1]

#write.table(seasonal_rain[, 'mosquitoes'], 'mosquito_seasonality.out', row.names=F, col.names=F, quote=F)

#max_rain = max(seasonal_rain$precip)

pdf('rain_refit2.pdf', width=20, height=6)
plot(sr_copy$daily, type='l', lwd=0.5)#, xlim=c(4000,5000))
lines(smooth.cropped, col=2, lwd=2)
#lines(seasonal_rain$smooth,col=2, lwd=2)
#lines(1:365+8, seasonal_rain$smooth,col=3, lwd=2)
dev.off()
