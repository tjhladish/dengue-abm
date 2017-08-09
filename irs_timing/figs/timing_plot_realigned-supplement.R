rm(list=ls())

require("RSQLite")
drv = dbDriver("SQLite")
db = dbConnect(drv, "./irs_timing-refit0.sqlite", flags=SQLITE_RO)

d <- dbGetQuery(db, 'select vector_control, timing, vc_coverage, campaign_duration, M.*
                      from par P, met M, job J
                      where P.serial = M.serial 
                    and P.serial = J.serial
                    and status = \'D\';')

npars = 4
serial_col = npars + 1
data_burnin = 5 # used 6 for timing plot
plot_years = 10  # used 5 for timing plot
last_col = serial_col + data_burnin + plot_years
month_starts = c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)
month_labels = c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

tags = d[,1:npars]
data = d[,(serial_col + data_burnin + 1):last_col]
#medians = aggregate(d[,5:34], by=list(d$vector_control, d$vc_coverage, d$timing), FUN=median)



db2 = dbConnect(drv, "rzero-refit.sqlite", flags=SQLITE_RO)
all_r0 = dbGetQuery(db2, 'select realization, m.* from met m, par p where m.serial=p.serial')
r0_100 = aggregate(all_r0, by=list(all_r0$realization), FUN=mean)
r0_100 = r0_100[,-c(1:3)]
r0 = sapply(r0_100, mean)
r0_q = sapply(r0_100, quantile, probs=c(0.025, 0.975))
r0_q025 = t(r0_q)[,1]
r0_q975 = t(r0_q)[,2]
smooth_tmp = smooth.spline(rep(r0, 3),spar=0.5)
r0_smooth = smooth_tmp$y[1:365 + 365]



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
mod3 = smooth.spline(rep(seasonal_rain$precip, 3),spar=0.6)
seasonal_rain$smooth = mod3$y[1:365 + 365]
seasonal_rain$mosquitoes = seasonal_rain$smooth/max(seasonal_rain$smooth)
# lag by one week to account for egg-to-emergence time, based on Table 2 in Rueda et al, Temperature-dependent
# development and survival rates of Culex quinquefasciatus and Aedes aegypti (Diptera: Culicidae), J Med Ent, 1990, 27:5
seasonal_rain$mosquitoes = seasonal_rain$mosquitoes[(0:364-7)%%365 + 1]

max_rain = max(seasonal_rain$precip)



















timing_effectiveness = function(tags, data, window, func, end) {
    if (end) {
        d = cbind(tags, rowSums(data[,(dim(data)[2]- window + 1):dim(data)[2]]))
    } else {
        d = cbind(tags, rowSums(data[,1:window]))
    }
    means = aggregate(d[,dim(d)[2]], by=list(d$vector_control, d$campaign_duration, d$vc_coverage, d$timing), FUN=func)
    names(means)=c('vc','dur','cov','day','sum')
    means$eff = 1 - means$sum/means$sum[1]
    return(means)
}

smoothed = function(y) {
    #x90 = (test$day[test$vc==1 & test$dur==dur_ & test$cov==cov_ & test$day < 365])
    #y90 = (test$eff[test$vc==1 & test$dur==dur_ & test$cov==cov_ & test$day < 365])
    #lines(x90, y90, lwd=0.5)
    y_3 <- c(y,y,y)
    f5 <- rep(1/5,5)
    y_mva_3 <- filter(y_3, f5, sides=2)
    y_mva = y_mva_3[54:106]
    #lines(x90, y_mva, lwd=lwd_, lty=lty_)
    return(y_mva)
}

rewrap_ts = function(x,y,delta){
    .tail = which(x+delta>365)
    x = c((x+delta)[.tail]-365, (x+delta)[-.tail])
    y = c(y[.tail], y[-.tail])
    return(list(x=x,y=y))
}

png('irs_timing-realigned-SI.png', width=2000, height=2200, res=300)
par(mfrow=c(3,1), mar=c(2.1,2.1,1,1), oma=c(3,3,0,0),las=1)
#dur_ = 1
end = F
for (cov_ in c(0.75, 0.5, 0.25)){
    test = timing_effectiveness(tags, data, plot_years, median, end)
    minval = min(c(test$eff[test$vc==1 & test$dur==0 & test$cov==cov_],test$eff[test$vc==1 & test$dur==1 & test$cov==cov_ & test$day < 365]))
    maxval = max(c(test$eff[test$vc==1 & test$dur==0 & test$cov==cov_],test$eff[test$vc==1 & test$dur==1 & test$cov==cov_ & test$day < 365]))
    .range = maxval-minval
    pad    = 0.05
    plot(0:365, ylim=c(minval-pad*.range,maxval+pad*.range), xlab='', ylab='', type='n', xaxt='n')
    axis(1,at=month_starts,labels = month_labels)
    
    lines((.range*r0_smooth/max(r0_smooth))+minval, col='#ff7700', lwd=2) # R0
    lines(1:365, (.range*seasonal_rain$mosquitoes/max(seasonal_rain$mosquitoes))+minval,col=4, lwd=1.5)          # mosquitoes
    
    .set = test$vc==1 & test$dur==0 & test$cov==cov_
    shifted = rewrap_ts(test$day[.set], test$eff[.set], 45)
    lines(shifted$x, shifted$y, type='l', lty=1, lwd=0.5)
    lines(shifted$x, smoothed(shifted$y), type='l', lty=1, lwd=2)

    .set = test$vc==1 & test$dur==1 & test$cov==cov_ & test$day < 365
    shifted = rewrap_ts(test$day[.set], test$eff[.set], 90)
    lines(shifted$x, shifted$y, type='l', lty=2, lwd=0.5)
    lines(shifted$x, smoothed(shifted$y), type='l', lty=2, lwd=2)
    # x90 = (test$day[])
    # y90 = (test$eff[test$vc==1 & test$dur==1 & test$cov==cov_ & test$day < 365])
    #lines(x90, y90, lty=2)
    abline(h=test$eff[test$vc==1 & test$dur==2 & test$cov==cov_], lty=3, lwd=2)
    #legend('bottomright',inset = c(0.08,0.16), legend = paste0(cov_*100, '% household coverage'), bty='n')
    text(160, minval, pos=4, labels = paste0(cov_*100, '% household coverage'), font=2)
    legend('topleft',title='Rollout period', legend = c('1 day', '90 days', '365 days'), lty=1:3, bty='n', lwd=2)
}
par(las=0)
mtext("Timing of peak coverage (Month)", side=1, outer=T, line=1)
mtext("Effectiveness (prevented cases)", side=2, outer=T, line=1)
#mtext("Effect of peak insecticide coverage on IRS effectiveness", side=3, outer=T)
dev.off()
