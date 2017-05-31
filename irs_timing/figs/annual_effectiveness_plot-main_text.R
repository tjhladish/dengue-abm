rm(list=ls())

require("RSQLite")
drv = dbDriver("SQLite")
# db = dbConnect(drv, "./irs_timing-refit0.sqlite", flags=SQLITE_RO)
#
# d <- dbGetQuery(db, 'select vector_control, timing, vc_coverage, campaign_duration, M.*
#                       from par P, met M, job J
#                       where P.serial = M.serial
#                     and P.serial = J.serial
#                     and status = \'D\';')

db = dbConnect(drv, "./irs_stopping-effect.sqlite", flags=SQLITE_RO)

d <- dbGetQuery(db, 'select vector_control, timing, vc_coverage, campaign_duration, strat_years, M.*
                      from par P, met M, job J
                      where P.serial = M.serial
                    and P.serial = J.serial
                    and status = \'D\';')



pop_size = 18.2 # in 100 thousands
npars = 5
serial_col = npars + 1
data_burnin = 5 # used 6 for timing plot
plot_years = 20  # used 5 for timing plot
last_col = serial_col + data_burnin + plot_years

tags = d[,1:npars]
#tags$strat_years = 50
data = d[,(serial_col + data_burnin + 1):last_col]

#tags_stop = d_stop[,1:npars]
#tags_stop$strat_years = 10
#data_stop = d_stop[,(serial_col + data_burnin + 1):last_col]
#tags_stop$timing = 154 ## DELETE THIS AFTER RE-RUNNING OTHER RESULTS
#medians = aggregate(d[,5:34], by=list(d$vector_control, d$vc_coverage, d$timing), FUN=median)

#npars = 5
#tags = rbind(tags, tags_stop)
#data = rbind(data, data_stop)

plot_effectiveness_over_time = function(tags, data, timing, plotcases=F, plotcumulative=F, plotcasesaverted=F) {
    medians = aggregate(data, by=list(tags$campaign_duration, tags$vc_coverage, tags$timing, tags$strat_years, tags$vector_control), FUN=median)
    names(medians)[1:npars]=c('dur','cov','day','strat','vc')
    medians = medians[medians$day==timing | medians$vc==0,]
    eff_vals = 1 - sweep(as.matrix(medians[-1,(npars+1):(dim(medians)[2])]), 2, as.numeric(medians[1,(npars+1):(dim(medians)[2])]), '/')
    eff = cbind(medians[-1,1:npars], eff_vals)
    #browser()
    cum_medians = cbind(medians[,1:npars], t(apply(medians[,-c(1:npars)], 1, cumsum)))
    cum_eff_vals = 1 - sweep(as.matrix(cum_medians[-1,(npars+1):(dim(cum_medians)[2])]), 2, as.numeric(cum_medians[1,(npars+1):(dim(cum_medians)[2])]), '/')
    cum_eff = cbind(cum_medians[-1,1:npars], cum_eff_vals)
    reds=c('#FF0000','#DD0000','#AA0000')
    .lwd = 1:3
    if (plotcases) {
        #browser()
        plotdata = rbind(medians[1,(npars+1):(npars+plot_years)], medians[medians$dur==1, (npars+1):(npars+plot_years)])/pop_size
        matplot(t(plotdata),
                lty=c(1,3:1), type='l', lwd=.lwd, col=c('black',reds),
                xlab='Year', ylab='Cases per 100,000 people', ylim=c(0,850), main='')
        legend('bottomright', legend=c('baseline','75% coverage','50% coverage','25% coverage'), lwd=.lwd, lty=c(1,1:3), bty='n', col=c('black',rev(reds)))
    } else if (plotcumulative) {
        matplot(t(cum_eff[cum_eff$dur==1, (npars+1):(npars+plot_years)]),
                lty=3:1, type='l', ylim=c(0,1), ylab='Cumulative effectiveness', lwd=.lwd, main='', xlab='Year', col=reds)
        legend('topright', legend=c('75% coverage','50% coverage','25% coverage'), lwd=.lwd, lty=1:3, bty='n', col=rev(reds))
    } else if (plotcasesaverted) {
        cum_ca_vals = -sweep(as.matrix(cum_medians[-1,(npars+1):(dim(cum_medians)[2])]), 2, as.numeric(cum_medians[1,(npars+1):(dim(cum_medians)[2])]), '-')
        cum_ca_vals = cum_ca_vals/(pop_size*1000)
        cum_ca = cbind(cum_medians[-1,1:npars], cum_ca_vals)
        matplot(t(cum_ca[cum_ca$dur==1, (npars+1):(npars+plot_years)]),
                lty=3:1, type='l', ylab='Cumulative cases averted per 100 people', lwd=.lwd, main='', xlab='Year', col=reds)
        legend('bottomright', legend=c('75% coverage','50% coverage','25% coverage'), lwd=.lwd, lty=1:3, bty='n', col=rev(reds))
    } else {
       # browser()
        matplot(t(eff[eff$dur==1 & eff$strat==10,(npars+1):(npars+plot_years)]),
                #lty=3:1, type='l', ylim=c(0,1.05), lwd=.lwd, ylab='Overall effectiveness', xlab='Year', col=reds)
                lty=3:1, type='l', lwd=.lwd, ylab='Overall effectiveness', xlab='Year', col='#999999')
        matlines(t(eff[eff$dur==1 & eff$strat==50, (npars+1):(npars+plot_years)]),
                lty=3:1, type='l', ylab='Cumulative cases averted per 100 people', lwd=.lwd, main='', xlab='Year', col=reds)
        abline(h=0,lty=2)
        legend('bottomright', legend=c('75% coverage','50% coverage','25% coverage'), lwd=rev(.lwd), lty=1:3, bty='n', col=rev(reds))
    }
    return(medians)
}

#png('vc_effectiveness.png', width=1200, height=1920, res=180)
png('vc_effectiveness-main-refit-20year-b.png', width=1500, height=1290, res=200)
#par(mfrow=c(2,2), mar=c(4.6,4.1,1.2,1), oma=c(0,0,2,0))

par(las=1,bty='L')
#plot_effectiveness_over_time(tags,data,182, plotcases=T)
plot_effectiveness_over_time(tags,data,152)
mtext('Simulated impact of IRS (90-day campaign, early June) over first 10 years',side = 3, line=2)
#plot_effectiveness_over_time(tags,data,182, plotcasesaverted=T)
#plot_effectiveness_over_time(tags,data,182, plotcumulative=T)
dev.off()
