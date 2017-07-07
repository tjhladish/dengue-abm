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

png('irs_timing-first10years-main.png', width=2500, height=860, res=225)
par(mfrow=c(1,3), mar=c(2.1,2.1,1,1), oma=c(3,3,3,0),las=1)
#dur_ = 1
end = F
for (cov_ in c(0.25, 0.5, 0.75)){
    test = timing_effectiveness(tags, data, plot_years, median, end)
    plot(0:365, ylim=c(0,1), xlab='', ylab='', type='n', xaxt='n')
    axis(1,at=month_starts,labels = month_labels)
    minval = min(c(test$eff[test$vc==1 & test$dur==0 & test$cov==cov_],test$eff[test$vc==1 & test$dur==1 & test$cov==cov_ & test$day < 365]))
    maxval = max(c(test$eff[test$vc==1 & test$dur==0 & test$cov==cov_],test$eff[test$vc==1 & test$dur==1 & test$cov==cov_ & test$day < 365]))
    rect(0,minval,365,maxval,border=NA,col='#cccccc')
    lines(test$day[test$vc==1 & test$dur==0 & test$cov==cov_], test$eff[test$vc==1 & test$dur==0 & test$cov==cov_],
          ylim=c(0,1), type='l', xlab='', ylab='', lty=1)
    
#    x90 = ((45+test$day[test$vc==1 & test$dur==1 & test$cov==cov_ & test$day < 365])%%365)[c(47:53,1:46)]
    x90 = (test$day[test$vc==1 & test$dur==1 & test$cov==cov_ & test$day < 365])
    y90 = (test$eff[test$vc==1 & test$dur==1 & test$cov==cov_ & test$day < 365])
    lines(x90, y90, lty=2)
    abline(h=test$eff[test$vc==1 & test$dur==2 & test$cov==cov_], lty=3)
    #legend('bottomright',inset = c(0.08,0.16), legend = paste0(cov_*100, '% household coverage'), bty='n')
    text(160, 0.185, pos=4, labels = paste0(cov_*100, '% household coverage'), font=2)
    legend('bottomright',legend = c('all houses treated in 1 day', 'houses treated across 90 days', 'houses treated continuously'), lty=1:3, bty='n')
}
par(las=0)
mtext("Campaign start (Month)", side=1, outer=T, line=1)
mtext("Effectiveness (prevented cases)", side=2, outer=T, line=1)
mtext("Effect of campaign start date on IRS effectiveness", side=3, outer=T)
dev.off()
