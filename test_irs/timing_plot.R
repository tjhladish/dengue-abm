rm(list=ls())

require("RSQLite")
drv = dbDriver("SQLite")
db = dbConnect(drv, "./irs_timing-1_vs_90_day_campaign.sqlite")

d <- dbGetQuery(db, 'select vector_control, timing, vc_coverage, campaign_duration, M.*
                      from parameters P, metrics M, jobs J
                      where P.serial = M.serial 
                    and P.serial = J.serial
                    and status = \'D\';')

npars = 4
serial_col = npars + 1
data_burnin = 0
plot_years = 30
last_col = serial_col + data_burnin + plot_years

tags = d[,1:npars]
data = d[,(serial_col + data_burnin + 1):last_col]
#medians = aggregate(d[,5:34], by=list(d$vector_control, d$vc_coverage, d$timing), FUN=median)

plot_effectiveness_over_time = function(tags, data, timing) {
    
    means = aggregate(data, by=list(tags$vector_control, tags$campaign_duration, tags$vc_coverage, tags$timing), FUN=mean)
    names(means)[1:npars]=c('vc','dur','cov','day')
    means = means[means$day==timing | means$vc==0,]
    eff_vals = 1 - sweep(as.matrix(means[-1,(npars+1):(dim(means)[2])]), 2, as.numeric(means[1,(npars+1):(dim(means)[2])]), '/')
    eff = cbind(means[-1,1:npars], eff_vals)
    
    # these aren't right:
    cum_means = cbind(means[,1:npars], t(apply(means[,-c(1:npars)], 1, cumsum)))
    cum_eff_vals = 1 - sweep(as.matrix(cum_means[-1,(npars+1):(dim(cum_means)[2])]), 2, as.numeric(cum_means[1,(npars+1):(dim(cum_means)[2])]), '/')
    cum_eff = cbind(cum_means[-1,1:npars], cum_eff_vals)
    
    par(mfrow=c(2,1), mar=c(2.1,2.1,1,1), oma=c(3,3,3,0))
    matplot(t(eff[eff$dur==90,5:29]), lty=1, type='l', ylim=c(-0.2, 0.45), ylab='Effectiveness')
    matplot(t(cum_eff[cum_eff$dur==90,5:29]), lty=1, type='l', ylim=c(-0.2, 0.45), ylab='Cumulative effectiveness')
    
    return(means)
}

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

png('ce_vc-equillibrium.png', width=6, height=8, res=150, units='in')
par(mfrow=c(3,1), mar=c(2.1,2.1,1,1), oma=c(3,3,3,0))
#dur_ = 1
end = F
for (cov_ in c(0.25, 0.5, 0.75)){
    test = timing_effectiveness(tags, data, 25, mean, end)
    plot(test$day[test$vc==1 & test$dur==1 & test$cov==cov_], test$eff[test$vc==1 & test$dur==1 & test$cov==cov_], ylim=c(0,0.5), type='l', xlab='', ylab='', lty=2)
    x90 = ((45+test$day[test$vc==1 & test$dur==90 & test$cov==cov_ & test$day < 365])%%365)[c(65:73,1:64)]
    y90 = (test$eff[test$vc==1 & test$dur==90 & test$cov==cov_ & test$day < 365])[c(65:73,1:64)]
    lines(x90, y90, ylim=c(0,0.5), type='l', xlab='', ylab='')
    legend('topleft',legend = paste0(cov_*100, '% household coverage'), bty='n')
    legend('topright',legend = c('all houses treated in 1 day', 'houses treated across 90 days'), lty=c(2,1), bty='n')
}
mtext("Campaign start (Julian day)", side=1, outer=T, line=1)
mtext("Effectiveness (prevented infections)", side=2, outer=T, line=1)
mtext("Effect of campaign timing on IRS effectiveness", side=3, outer=T)
dev.off()
