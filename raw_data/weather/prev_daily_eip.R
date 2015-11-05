all_data = read.csv("NOAA_yucatan_daily.csv")

# work with Merida data only
d = all_data[all_data$NAME=='AEROP.INTERNACIONAL',]
d$DATE = as.Date(d$DATE)
# temps currently in tenths of degrees
d$TMAX = d$TMAX/10
d$TMIN = d$TMIN/10

drel <- data.table(d)[,list(TMAX, TMIN),keyby=DATE]

alpha = 0.001
validLows = quantile(d$TMIN,  probs = c(alpha, 1-alpha), na.rm=T)
validHighs = quantile(d$TMAX,  probs = c(alpha, 1-alpha), na.rm=T)
outlier_plots = function(save=F) {
  if(save) png("outlier_temps.png", width=2000, height=1600, res=200)
  oldpar <- par(mfrow=c(2,2))
  
  tmp_y_data = head(sort(na.omit(d$TMIN)), n=100)
  plot(tmp_y_data, ylab='Temp (deg C)', main='Smallest daily low temperatures', cex=0.5, 
       col=1+as.numeric(tmp_y_data < validLows[1]))
  
  tmp_y_data = head(sort(na.omit(d$TMIN), decreasing = T), n=100)
  plot(tmp_y_data, ylab='Temp (deg C)', main='Largest daily low temperatures', cex=0.5, 
       col=1+as.numeric(tmp_y_data > validLows[2]))
  
  tmp_y_data = head(sort(na.omit(d$TMAX)), n=100)
  plot(tmp_y_data, ylab='Temp (deg C)', main='Smallest daily high temperatures', cex=0.5, 
       col=1+as.numeric(tmp_y_data < validHighs[1]))
  
  tmp_y_data = head(sort(na.omit(d$TMAX), decreasing = T), n=100)
  plot(tmp_y_data, ylab='Temp (deg C)', main='Largest daily high temperatures', cex=0.5, 
       col=1+as.numeric(tmp_y_data > validHighs[2]))
  
  if (save) dev.off()
  par(oldpar)
}

outlier_plots()

# throw out most extreme 0.1%
d$TMIN[d$TMIN < validLows[1]] = NA
d$TMIN[d$TMIN > validLows[2]] = NA
d$TMAX[d$TMAX < validHighs[1]] = NA
d$TMAX[d$TMAX > validHighs[2]] = NA

# generate reconstructed temps data structure, starting with date sequence
# TMINr and TMAXr will be the reconstructed temps
rtemps=data.table(DATE=seq(as.Date("1979/1/1"), as.Date("2014/1/1")-1, "days"), TMINr=NA, TMAXr=NA)
rtemps$TMIN = d$TMIN[match(rtemps$DATE, d$DATE)]
rtemps$TMAX = d$TMAX[match(rtemps$DATE, d$DATE)]

# grab miami data in order to reconstruct merida data
m = read.table("NOAA_miami.csv", header=T)
m$DATE = as.Date(paste(m$YEAR,"/",m$MONTH, "/", m$DAY, sep=""))
rtemps$TMINmia = m$TMIN[match(rtemps$DATE, m$DATE)]/10
rtemps$TMAXmia = m$TMAX[match(rtemps$DATE, m$DATE)]/10
rtemps$month = as.factor(format(rtemps$DATE, '%b'))

# regressions using miami temp and month
reg_min = lm(TMIN ~ TMINmia + month, data=rtemps)
reg_max = lm(TMAX ~ TMAXmia + month, data=rtemps)

# predict missing values more sensibly
rtemps$TMINr[is.na(rtemps$TMIN)] = predict(reg_min, newdata=rtemps[is.na(rtemps$TMIN),], type='response')
rtemps$TMAXr[is.na(rtemps$TMAX)] = predict(reg_max, newdata=rtemps[is.na(rtemps$TMAX),], type='response')

# store fahrenheits, too, for convenience for us non-metric weirdos
#d$TMINf = d$TMIN * 9/5 + 32
#d$TMAXf = d$TMAX * 9/5 + 32

# daily mean calculated as simple average of daily min and max
#d$TMEAN = (d$TMIN + d$TMAX)/2
#d$TMEANf = (d$TMINf + d$TMAXf)/2

# generate date sequence
# pull out means for the four locations we have, can't remember what the match part does exactly
# We're only really using the merida values currently

min_max_plot = function(save=F) {
  if(save) pdf("Merida_min_max_temps.pdf", width=11, height=8.5)
  #png("reconstructed_merida_min_max_temperatures.png", width=1800, height=800, res=200)
  oldpar <- par(mfrow=c(2,1), mar=c(2.5,4,3,1))
  # merida mins + predictions from miami mins
  plot(rtemps$DATE, rtemps$TMIN, type='l', col='black', ylab="", xlab="", main="Mérida daily low temperatures")
  points(rtemps$DATE, rtemps$TMINr, type='l', col='blue')
  # merida maxes + predictions from miami maxes
  plot(rtemps$DATE, rtemps$TMAX, type='l', col='black', ylab="", xlab="", main="Mérida daily high temperatures")
  points(rtemps$DATE, rtemps$TMAXr, type='l', col='red')
  if (save) dev.off()
  par(oldpar)
}

min_max_plot()

# create new columns that merge the merida data and the reconstructed values
rtemps$TMINmerged = rtemps$TMIN
rtemps$TMINmerged[is.na(rtemps$TMINmerged)] = rtemps$TMINr[is.na(rtemps$TMINmerged)]
rtemps$TMAXmerged = rtemps$TMAX
rtemps$TMAXmerged[is.na(rtemps$TMAXmerged)] = rtemps$TMAXr[is.na(rtemps$TMAXmerged)]

# EIP log-normal function from http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0050972
eip = function(T) { exp(exp(2.9 -0.08*T) + 1/(2*4.9)) }

# EIR is extrinsic incubation rate, 1/EIP
# assumes square temp function, i.e. 12 hours at daily max temp and 12 hours at daily min temp
rtemps$EIR = 1/rowMeans(data.frame(eip(rtemps$TMINmerged), eip(rtemps$TMAXmerged)))

smoothed_eip = function(current_day) {
  total = 0
  eip = 0
  while (total < 1) {
    total = total + rtemps$EIR[current_day]
    eip = eip + 1
    if (current_day < length(rtemps$EIR)) {current_day = current_day + 1}
  }
  return(eip)
}

for (day in 1:dim(rtemps)[1]) { 
  rtemps[day, 'EIP'] = smoothed_eip(day) 
}

eip_plots = function() {
  pdf("Merida_EIP_ts.pdf", width=11, height=8.5)
  #plot(rtemps$DATE, 1/rtemps$EIR, log='y', type='l', xlab='Date', ylab='EIP', main='Estimated daily EIP for Mérida')
  plot(rtemps$DATE, rtemps$EIP, type='l', xlab='Date', ylab='EIP', main='Estimated daily EIP for Mérida')
  dev.off()

  pdf("Merida_EIP_hist.pdf", width=11, height=8.5)
  hist(rtemps$EIP, breaks=100, axes=F, xlab='EIP', main='Histogram of daily EIP values')
  axis(2)
  axis(1,at=seq(10,90,10))
  dev.off()
}

#min_max_plot()
#eip_plots()
