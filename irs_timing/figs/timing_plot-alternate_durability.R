rm(list=ls())

require("RSQLite")
drv = dbDriver("SQLite")
db1 = dbConnect(drv, "./irs_timing-refit0.sqlite", flags=SQLITE_RO)
db2 = dbConnect(drv, "./irs_insecticide-durability-effect.sqlite", flags=SQLITE_RO)

d1 <- dbGetQuery(db1, 'select vector_control, timing, vc_coverage, 90 as eff_days, M.*
                      from par P, met M, job J
                      where P.serial = M.serial 
                    and campaign_duration = 1
                    and vc_coverage = 0.75
                    and P.serial = J.serial
                    and status = \'D\';')

d2 <- dbGetQuery(db2, 'select vector_control, timing, vc_coverage, eff_days, M.*
                      from par P, met M, job J
                      where P.serial = M.serial 
                    and P.serial = J.serial
                    and status = \'D\';')
d1 = d1[,1:20]
d2 = d2[,1:20]

d = rbind(d1,d2)

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
    means = aggregate(d[,dim(d)[2]], by=list(d$vector_control, d$vc_coverage, d$timing, d$eff_days), FUN=func)
    names(means)=c('vc','cov','day','eff_days','sum')
    means$eff = 1 - means$sum/means$sum[1]
    return(means)
}

png('irs_timing-insecticide_durability.png', width=1500, height=1000, res=200)
par(mfrow=c(1,1))
par(mar=c(5.1, 4.1, 1.1, 2.1))
#par(oma=c(2.5,3,0.5,0.5))
par(las=1)

#par(mfrow=c(2,1), mar=c(2.1,2.1,1,1), oma=c(3,3,3,0),las=1)
end = F


x_lab_idx = seq(1,12,2)
test = timing_effectiveness(tags, data, plot_years, median, end)

plot(0:365, ylim=c(0,1), xlab="Campaign start (Month)", ylab="Effectiveness (prevented cases)", type='n', xaxt='n', yaxt='n')
axis(1,at=month_starts[x_lab_idx],labels = month_labels[x_lab_idx])
axis(1,at=month_starts[-x_lab_idx], labels = F)
axis(2,at=seq(0,1,0.2))
axis(2,at=seq(0.1,0.9,0.2),labels=F)

plot_effectiveness_curve = function(cov_, lwd_, eff_days_) {
    lty_=1
    #if (dur_==0) lty_=2
    x90 = (test$day[test$vc==1 & test$eff_days==eff_days_ & test$cov==cov_ & test$day < 365])
    y90 = (test$eff[test$vc==1 & test$eff_days==eff_days_ & test$cov==cov_ & test$day < 365])
    lines(x90, y90, lwd=0.5)
    y90_3 <- c(y90,y90,y90)
    f5 <- rep(1/5,5)
    y_mva_3 <- filter(y90_3, f5, sides=2)
    y_mva = y_mva_3[54:106]
    lines(x90,y_mva, lwd=lwd_, lty=lty_)
    print(c(x90[which.max(y_mva)],max(y_mva), x90[which.min(y_mva)], min(y_mva)))
}

#text(0, 0.03, pos=4, labels = 'Houses treated across 90 days', font=2)
text(0, 0.97, pos=4, labels = '90 day campaigns, 75% coverage', font=2)
legend('topright',legend = c('150-day durability', '90-day durability', '30-day durability'), lwd=4:2, bty='n')


cov_ = 0.75
ylim_ = c(0,1)
#plot(0:365, ylim=ylim_, xlab='', ylab='', type='n', xaxt='n')
axis(1,at=month_starts[x_lab_idx],labels = month_labels[x_lab_idx])
axis(1,at=month_starts, labels = F)

plot_effectiveness_curve(cov_ = cov_, lwd_ = 2, eff_days_ = 30)
plot_effectiveness_curve(cov_ = cov_, lwd_ = 3, eff_days_ = 90)
plot_effectiveness_curve(cov_ = cov_, lwd_ = 4, eff_days_ = 150)

#text(0, ylim_[2]-(0.05*(ylim_[2]-ylim_[1])), pos=4, labels = paste0(cov_*100, '% household coverage'), font=2)
#legend('topright',legend = c('all houses treated in 1 day', 'houses treated across 90 days', 'houses treated continuously'), lty=c(2,1,3), lwd=lwd_, bty='n')

par(las=0)

#mtext("Effect of campaign start date on IRS effectiveness", side=3, outer=T)
dev.off()

