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
month_starts = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366)
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

png('irs_timing-first10years-75-main-refit.png', width=1500, height=1500, res=200)
par(mfrow=c(2,1))
par(mar=c(2.0, 2.1, 0, 1.1))
par(oma=c(1.5,3,0.5,0.5))
par(las=1)

#par(mfrow=c(2,1), mar=c(2.1,2.1,1,1), oma=c(3,3,3,0),las=1)
end = F


x_lab_idx = seq(1,12,2)
test = timing_effectiveness(tags, data, plot_years, median, end)
plot(1:365, ylim=c(0,1), xlab='', ylab='', type='n', xaxt='n', yaxt='n')
axis(1,at=month_starts[x_lab_idx],labels = month_labels[x_lab_idx], lwd = 0, line=-0.3)
axis(1,at=month_starts, labels = F)
axis(2,at=seq(0,1,0.2))
axis(2,at=seq(0.1,0.9,0.2),labels=F)

plot_effectiveness_curve = function(cov_, lwd_, dur_=1) {
    lty_=1
    if (dur_==0) lty_=2
    x90 = (test$day[test$vc==1 & test$dur==dur_ & test$cov==cov_ & test$day < 365])
    y90 = (test$eff[test$vc==1 & test$dur==dur_ & test$cov==cov_ & test$day < 365])
    lines(x90, y90, lwd=0.5)
    y90_3 <- c(y90,y90,y90)
    f5 <- rep(1/5,5)
    y_mva_3 <- filter(y90_3, f5, sides=2)
    y_mva = y_mva_3[54:106]
    lines(x90,y_mva, lwd=lwd_, lty=lty_)
    print(c(x90[which.max(y_mva)],max(y_mva), x90[which.min(y_mva)], min(y_mva)))
}

plot_effectiveness_curve(0.25,1)
plot_effectiveness_curve(0.5,2)
plot_effectiveness_curve(0.75,3)

#text(0, 0.03, pos=4, labels = 'Houses treated across 90 days', font=2)
text(0.04*365, 0.97, pos=4, labels = 'Houses treated across 90 days', font=2)
text(0, 0.97, cex=2, font=2,labels='c')
legend('topright',legend = c('75% coverage', '50% coverage', '25% coverage'), lwd=3:1, bty='n')


cov_ = 0.75
ylim_ = c(0.4,1)
lwd_ = 3 #if (cov_==0.5) 2 else 3
plot(1:365, ylim=ylim_, xlab='', ylab='', type='n', xaxt='n')
axis(1,at=month_starts[x_lab_idx],labels = month_labels[x_lab_idx], lwd = 0, line=-0.3)
axis(1, at=month_starts, labels = F)
minval = min(c(test$eff[test$vc==1 & test$dur==0 & test$cov==cov_],test$eff[test$vc==1 & test$dur==1 & test$cov==cov_ & test$day < 365]))
maxval = max(c(test$eff[test$vc==1 & test$dur==0 & test$cov==cov_],test$eff[test$vc==1 & test$dur==1 & test$cov==cov_ & test$day < 365]))
#rect(0,minval,365,maxval,border=NA,col='#cccccc')

plot_effectiveness_curve(cov_ = cov_, lwd_ = lwd_, dur_ = 0)
plot_effectiveness_curve(cov_ = cov_, lwd_ = lwd_)
continuous_Tx_eff = test$eff[test$vc==1 & test$dur==2 & test$cov==cov_]
lines(c(1,365), c(continuous_Tx_eff, continuous_Tx_eff), lty=3, lwd=lwd_)
#abline(h=, lty=3, lwd=lwd_)
text(0.04*365, ylim_[2]-(0.03*(ylim_[2]-ylim_[1])), pos=4, labels = paste0(cov_*100, '% household coverage'), font=2)
text(0, ylim_[2]-(0.03*(ylim_[2]-ylim_[1])), cex=2, font=2,labels='d')
#text(0, ylim_[1]+(0.03*(ylim_[2]-ylim_[1])), pos=4, labels = paste0(cov_*100, '% household coverage'), font=2)
legend('topright',legend = c('all houses treated in 1 day', 'houses treated across 90 days', 'houses treated across 365 days'), lty=c(2,1,3), lwd=lwd_, bty='n')

par(las=0)

mtext("Campaign start (Month)", side=1, outer=T, line=0)
mtext("Effectiveness (prevented cases)", side=2, outer=T, line=1)
#mtext("Effect of campaign start date on IRS effectiveness", side=3, outer=T)
dev.off()
