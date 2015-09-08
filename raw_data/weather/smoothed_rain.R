all_data = read.csv("NOAA_yucatan_daily.csv")

# work with Merida data only
d = all_data[all_data$NAME=='AEROP.INTERNACIONAL',]
d$DATE = as.Date(d$DATE)

precip=data.frame(DATE=seq(as.Date("1979/1/1"), as.Date("2014/1/1"), "days"))
precip$PRCP = d$PRCP[match(precip$DATE, d$DATE)]

precip$bool0 = precip$PRCP > 0
precip$month = strftime(precip$DATE, '%m')
precip$day = strftime(precip$DATE, '%m-%d')
seasonal_rain = aggregate(precip$bool0, by=list(precip$day), mean, na.rm=T)

names(seasonal_rain) = c('day', 'precip')
seasonal_rain = seasonal_rain[-which(seasonal_rain$day == '02-29'),]

julian_days = 1:365
#lo = loess(seasonal_rain$precip~julian_days)
#plot(seasonal_rain$precip)
#lines(predict(lo), col='blue', lwd=2)

#x = ts(seasonal_rain$precip)
#fit = ts(tsSmooth(StructTS(x))[,-2])
#tsp(fit) = tsp(x)

#pdf('rain_fit.pdf', width=10, height=6)
#plot(x)
#lines(fit, col=2)
#dev.off()

#forward=fft(seasonal_rain$precip); forward[ abs(forward)<5 ]=0;  smoothed=fft(forward, inverse=TRUE);
#pdf('rain_fit.pdf', width=10, height=6); plot(x); lines(Re(smoothed/length(smoothed)),col=2); dev.off()

mod3 = smooth.spline(rep(seasonal_rain$precip, 3),spar=0.6)
seasonal_rain$smooth = mod3$y[1:365 + 365]
pdf('rain_fit.pdf', width=10, height=6); plot(seasonal_rain$precip, type='l', lwd=0.5); lines(seasonal_rain$smooth,col=2, lwd=2); lines(1:365+8, seasonal_rain$smooth,col=3, lwd=2); dev.off()

write.table(seasonal_rain[, c('smooth', 'precip', 'day')], 'merida_preciptiation.out', row.names=F, col.names=F, quote=F) 

seasonal_rain$mosquitoes = seasonal_rain$smooth/max(seasonal_rain$smooth)
# lag by one week to account for egg-to-emergence time, based on Table 2 in Rueda et al, Temperature-dependent
# development and survival rates of Culex quinquefasciatus and Aedes aegypti (Diptera: Culicidae), J Med Ent, 1990, 27:5
seasonal_rain$mosquitoes = seasonal_rain$mosquitoes[(0:364-7)%%365 + 1]

write.table(seasonal_rain[, 'mosquitoes'], 'mosquito_seasonality.out', row.names=F, col.names=F, quote=F)

#(1:10+3)%%11
max_rain = max(seasonal_rain$precip)
e = read.table('seasonal_eip.out', header=T)
e_avg = read.table('seasonal_avg_eip.out', header=T)
rescaled_e_avg = e_avg[,1]/max(e_avg[,1])

t = read.table('merida_temps.out', header=F)
rescaled_tmin = e[,3]/max(t[,1])
rescaled_tmax = e[,4]/max(t[,1])
rescaled_t    = t[,1]/max(t[,1])
#rain = c(0.189, 0.179,0.128,0.123,0.0956,0.195,0.777,0.940,0.901,1.0,0.491,0.301,0.199,0.189)
month_days = c(31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365) - 15
month_days2 = c(1,month_days,365)
rescaled_rain = rain/max(rain)
pdf('rain_fit.pdf', width=10, height=6); 
plot(seasonal_rain$precip/max(seasonal_rain$smooth), type='l', lwd=0.5, ylab='Scale', xlab='Julian day', main='Seasonality in Merida, Yucatan'); 
lines(seasonal_rain$smooth/max(seasonal_rain$smooth),col=4, lwd=2);# lines(1:365+8, seasonal_rain$smooth,col=3, lwd=2); dev.off()
lines(rescaled_t, col=2, lwd=2); 
lines(rescaled_tmin, col=2);
lines(rescaled_tmax, col=2); 
lines(rescaled_e_avg, col=3, lwd=2); 
points(month_days2, rescaled_rain, type='l', col='blue', lty=2); 
legend('topright', legend=c('Mean temp', 'Min/max temp', 'EIP', 'Pr{precip}', 'Smoothed precip', 'Old precip model'), col=c('red','red','green','black','blue','blue'), lwd=c(2,1,2,0.5,2,2), lty=c(1,1,1,1,1,2), bty='n', cex=0.8)
dev.off()
