require(data.table)
require(lubridate)

dateRange <- list(start=as.Date("1979/1/1"), end=as.Date("2014/1/1")-1)

miami  <- readRDS("miami.RData")[  between(DATE, dateRange$start, dateRange$end) ]
merida <- readRDS("merida.RData")[ between(DATE, dateRange$start, dateRange$end) ]

require(zoo)
require(tsoutliers)

tm.ts <- zoo(merida$TMAX, merida$DATE)

extremeLimits <- function(dt, alpha = 0.0005, ps=c(alpha, 1-alpha)) with(dt,{
  lows <- quantile(TMIN, probs=ps, na.rm=T)
  highs <- quantile(TMAX, probs=ps, na.rm=T)
  list(lows = lows, highs = highs)
})

valid = extremeLimits(merida)

require(ggplot2)
require(reshape2)
require(grid)

showValids <- function(dt, validation, ...) with(validation, {
  molten <- melt(dt, id.var="DATE")[!is.na(value)][, src := factor("VALID", levels=c("VALID","LOW","HIGH","REG."), ordered = T)][, year := year(DATE) ]
  molten[variable == "TMIN" & value < lows[1], src := factor("LOW", levels=c("VALID","LOW","HIGH","REG."), ordered = T)]
  molten[variable == "TMIN" & value > lows[2], src := factor("HIGH", levels=c("VALID","LOW","HIGH","REG."), ordered = T)]
  molten[variable == "TMAX" & value < highs[1], src := factor("LOW", levels=c("VALID","LOW","HIGH","REG."), ordered = T)]
  molten[variable == "TMAX" & value > highs[2], src := factor("HIGH", levels=c("VALID","LOW","HIGH","REG."), ordered = T)]
  ggplot(molten) + theme_bw() + theme(...) + aes(color=variable, x=yday(DATE), y=value) +
    geom_line(mapping = aes(group=year), alpha=0.1) + geom_point(mapping=aes(shape=src), data=molten[src!="VALID"]) + facet_grid(. ~ variable) +
    scale_color_manual(values=c(TMAX="red", TMIN="blue")) +
    scale_alpha_manual(values=c(VALID=0.2, LOW=0.5, HIGH=0.5, "REG."=1))
})

plt <- showValids(merida, valid, panel.border=element_blank(), panel.margin=unit(0.1, "lines"))

# throw out most extreme alpha*2
meridacensored = with(valid, merida[between(TMIN, lows[1], lows[2], incbounds = FALSE) & between(TMAX, highs[1], highs[2], incbounds = FALSE)])

reconstructedtemps = merge(
    data.table(DATE=seq(dateRange$start, dateRange$end, "days"), key="DATE"),
    meridacensored,
    by="DATE", all.x=T
)

# TODO switch to principle components regression?  cbind(TMAX, TMIN) ~ i.TMAX + i.TMIN + month (+ year?)
# regressions using miami temp and month
regressors = reconstructedtemps[miami][, month := as.factor(format(DATE, '%b'))]
reg = with(regressors,{
  list(
    min=lm(TMIN ~ i.TMIN + month),
    max=lm(TMAX ~ i.TMAX + month)
  )
})

reconstructedtemps[is.na(TMAX), TMAX := regressors[is.na(TMAX), predict(reg$max, newdata=.SD)] ]
reconstructedtemps[is.na(TMIN), TMIN := regressors[is.na(TMIN), predict(reg$min, newdata=.SD)] ]

showReconstructed <- function(regs, plt) {
  dt <- regs[is.na(TMAX) | is.na(TMIN), list(TMAX=i.TMAX, TMIN=i.TMIN, src="REG."), keyby=DATE]
  molten <- melt(dt, id.vars = c("DATE","src"))[, year := year(DATE)]
  plt + geom_point(data=molten)
}

showReconstructed(regressors, plt)

censoredtemps = reconstructedtemps[!(month(DATE)==2 & day(DATE)==29)]

## at this point, slight differences in reconstructed values, due to initial (here) vs later (in orig) pruning by date

# from Chan-Johansson
beta0 = 2.9
betaT = -0.08
var = 1/4.9

eip = function(tempC, beta0, betaT, var) { exp(exp(beta0 + betaT*tempC) + var/2) }

compareEIPs <- function(dt) {
  
}

eir = function(tempC, beta0, betaT, var) { 1/exp(exp(beta0 + betaT*tempC) + var/2) }

plot(
  censoredtemps[,
    list(eip=(eip(TMAX, beta0, betaT, var)+eip(TMIN, beta0, betaT, var))/2),
    keyby=DATE
  ][, mean(eip), keyby=as.Date(paste("2001", month(DATE), mday(DATE),sep="-"))],
  ylim=c(0,1000)
)

eirs = reconstructedtemps[,
  list(
    EIRMIN=eir(TMIN,beta0,betaT,var)/2,
    EIRMAX=eir(TMAX,beta0,betaT,var)/2
  ),                  
  keyby=DATE
]

# look backwards to find upper limit on start date
limdate = eirs[dim(eirs)[1]:1,list(datelim=cumsum(EIRMAX+EIRMIN), DATE)][datelim > 1][1,DATE]-1

lookaheads = eirs[DATE <= limdate, {
  cumeir = EIRMAX/2 + EIRMIN # assume initial bite is during the day
  ahead = DATE
  while(cumeir < 1) {
    ahead = ahead+1
    slice = eirs[DATE == ahead]
    add = slice$EIRMAX
    cumeir = cumeir + add
    if (cumeir < 1) {
      add = slice$EIRMIN
      cumeir = cumeir + add
    }
  }
  del = cumeir - 1
  as.numeric(ahead-DATE) + (1-del/add)*0.5
}, keyby=DATE]

eipmeantoTeff = function(meanEIP, beta0, betaT, var) (log(log(meanEIP) - var/2) - beta0)/betaT

effeips = lookaheads[!(month(DATE)==2 & mday(DATE)==29),
  mean(V1),
  keyby=list(month(DATE),mday(DATE))
][,list(eipeff=V1), keyby=as.Date(paste("2001",month,mday,sep="-"))]



# store fahrenheits, too, for convenience for us non-metric weirdos
#d$TMINf = d$TMIN * 9/5 + 32
#d$TMAXf = d$TMAX * 9/5 + 32

# daily mean calculated as simple average of daily min and max
#d$TMEAN = (d$TMIN + d$TMAX)/2
#d$TMEANf = (d$TMINf + d$TMAXf)/2

# generate date sequence
# pull out means for the four locations we have, can't remember what the match part does exactly
# We're only really using the merida values currently

min_max_plot = function() {
  pdf("Merida_min_max_temps.pdf", width=11, height=8.5)
  #png("reconstructed_merida_min_max_temperatures.png", width=1800, height=800, res=200)
  par(mfrow=c(2,1))
  par(mar=c(2.5,4,3,1))
  # merida mins + predictions from miami mins
  plot(rtemps$DATE, rtemps$TMIN, type='l', col='black', ylab="", xlab="", main="Mérida daily low temperatures")
  points(rtemps$DATE, rtemps$TMINr, type='l', col='blue')
  # merida maxes + predictions from miami maxes
  plot(rtemps$DATE, rtemps$TMAX, type='l', col='black', ylab="", xlab="", main="Mérida daily high temperatures")
  points(rtemps$DATE, rtemps$TMAXr, type='l', col='red')
  dev.off()
}

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

mean_seasonal_eip = aggregate(rtemps$EIP, by=list(as.factor(format(rtemps$DATE, '%m-%d'))), mean)
names(mean_seasonal_eip)=c('day', 'eip')

mean_seasonal_eip = mean_seasonal_eip[-which(mean_seasonal_eip$day=='02-29'),]
mean_seasonal_eip$eip = round(mean_seasonal_eip$eip)
write.table(mean_seasonal_eip[,c(2,1)], 'seasonal_avg_eip.out', quote=F, row.names=F, col.names=F)


# mean_seasonality = aggregate(rtemps[c('TMINmerged', 'TMAXmerged')], by=list(as.factor(format(rtemps$DATE, '%m-%d'))), mean)
# names(mean_seasonality) = c('day', 'TMIN', 'TMAX')
# mean_seasonality$EIR = 1/rowMeans(data.frame(eip(mean_seasonality$TMIN), eip(mean_seasonality$TMAX)))
#
# smoothed_eip_w_wrapping = function(current_day) {
#     total = 0
#     eip = 0
#     while (total < 1) {
#         total = total + mean_seasonality$EIR[current_day]
#         eip = eip + 1
#         if (current_day < length(mean_seasonality$EIR)) {
#             current_day = current_day + 1
#         } else {
#             current_day = 1
#         }
#     }
#     return(eip)
# }
#
# for (day in 1:dim(mean_seasonality)[1]) {
#     mean_seasonality[day, 'EIP'] = smoothed_eip_w_wrapping(day)
# }
#
# # compare averaging temps and then calculating EIP with calculating daily EIP and then averaging
# #plot(mean_seasonality$EIP, type='l', xlab='Julian day', ylab='EIP')
# #lines(round(mean_seasonal_eip$eip), type='l', col='blue')
# #legend('topleft', legend=c('EIP(average temps)', 'average(EIP(daily temps))'), fill=c('black','blue'), bty='n')
#
# # get rid of leap day, since sim expects 365 days/year
# mean_seasonality = mean_seasonality[-which(mean_seasonality$day=='02-29'),]
# write.table(mean_seasonality[,c(5,1,2,3,4)], 'seasonal_eip.out', quote=F, row.names=F)

outlier_plots = function() {
    png("outlier_temps.png", width=2000, height=1600, res=200)
    par(mfrow=c(2,2))
    
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
    
    dev.off()
}
