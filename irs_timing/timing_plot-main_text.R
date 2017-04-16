rm(list=ls())

require("RSQLite")
drv = dbDriver("SQLite")
db = dbConnect(drv, "./irs_timing-refit0.sqlite")

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

png('irs_timing-first10years-75-main-refit.png', width=1000, height=1200, res=150)
par(mfrow=c(2,1), mar=c(2.1,2.1,1,1), oma=c(3,3,3,0),las=1)
end = F


x_lab_idx = seq(1,12,2)
test = timing_effectiveness(tags, data, plot_years, median, end)
plot(0:365, ylim=c(0,1), xlab='', ylab='', type='n', xaxt='n', yaxt='n')
axis(1,at=month_starts[x_lab_idx],labels = month_labels[x_lab_idx])
axis(1,at=month_starts[-x_lab_idx], labels = F)
axis(2,at=seq(0,1,0.2))
axis(2,at=seq(0.1,0.9,0.2),labels=F)

plot_effectiveness_curve = function(cov_, lwd_) {
    x90 = (test$day[test$vc==1 & test$dur==1 & test$cov==cov_ & test$day < 365])
    y90 = (test$eff[test$vc==1 & test$dur==1 & test$cov==cov_ & test$day < 365])
    points(x90, y90, pch=20, cex=0.5)
    y90_3 <- c(y90,y90,y90)
    x.mid <- seq(0,364,7); offset <- 371
    y90.smooth <- lowess(y90_3, f=1/20)
    y90.smoothed_vals = y90.smooth$y[(x.mid+offset)/7]
    lines(x.mid, y90.smooth$y[(x.mid+offset)/7], lwd=lwd_)
    #print(c(min(y90.smoothed_vals),max(y90.smoothed_vals)))
    print(c(y90.smooth$y[22+53], y90.smooth$y[44+53]))
    #browser()
}

plot_effectiveness_curve(0.25,1)
plot_effectiveness_curve(0.5,2)
plot_effectiveness_curve(0.75,3)

#text(0, 0.03, pos=4, labels = 'Houses treated across 90 days', font=2)
text(0, 0.95, pos=4, labels = 'Houses treated across 90 days', font=2)
legend('topright',legend = c('75% coverage', '50% coverage', '25% coverage'), lwd=3:1, bty='n')


cov_ = 0.75
ylim_ = if (cov_==0.5) c(0.35,0.75) else c(0.4,1)
lwd_ = 2 #if (cov_==0.5) 2 else 3
plot(0:365, ylim=ylim_, xlab='', ylab='', type='n', xaxt='n')
axis(1,at=month_starts[x_lab_idx],labels = month_labels[x_lab_idx])
axis(1,at=month_starts, labels = F)
minval = min(c(test$eff[test$vc==1 & test$dur==0 & test$cov==cov_],test$eff[test$vc==1 & test$dur==1 & test$cov==cov_ & test$day < 365]))
maxval = max(c(test$eff[test$vc==1 & test$dur==0 & test$cov==cov_],test$eff[test$vc==1 & test$dur==1 & test$cov==cov_ & test$day < 365]))
rect(0,minval,365,maxval,border=NA,col='#cccccc')
lines(test$day[test$vc==1 & test$dur==0 & test$cov==cov_], test$eff[test$vc==1 & test$dur==0 & test$cov==cov_],
      ylim=c(0,1), type='l', xlab='', ylab='', lty=2, lwd=lwd_)

x90 = (test$day[test$vc==1 & test$dur==1 & test$cov==cov_ & test$day < 365])
y90 = (test$eff[test$vc==1 & test$dur==1 & test$cov==cov_ & test$day < 365])
lines(x90, y90, lty=1, lwd=lwd_)
abline(h=test$eff[test$vc==1 & test$dur==2 & test$cov==cov_], lty=3, lwd=lwd_)
text(0, ylim_[2]-(0.05*(ylim_[2]-ylim_[1])), pos=4, labels = paste0(cov_*100, '% household coverage'), font=2)
#text(0, ylim_[1]+(0.03*(ylim_[2]-ylim_[1])), pos=4, labels = paste0(cov_*100, '% household coverage'), font=2)
legend('topright',legend = c('all houses treated in 1 day', 'houses treated across 90 days', 'houses treated continuously'), lty=c(2,1,3), lwd=lwd_, bty='n')

par(las=0)

mtext("Campaign start (Month)", side=1, outer=T, line=1)
mtext("Effectiveness (prevented cases)", side=2, outer=T, line=1)
mtext("Effect of campaign start date on IRS effectiveness", side=3, outer=T)
dev.off()
