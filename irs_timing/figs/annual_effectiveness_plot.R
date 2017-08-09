rm(list=ls())

require("RSQLite")
drv = dbDriver("SQLite")
db = dbConnect(drv, "./irs_stopping-effect_rerun.sqlite", flags=SQLITE_RO)

d <- dbGetQuery(db, 'select vector_control, timing, vc_coverage, campaign_duration, M.*
                      from par P, met M, job J
                      where P.serial = M.serial
                    and strat_years = 50
                    and (vc_coverage = 0.75 or vector_control)
                    and P.serial = J.serial
                    and status = \'D\';')

#d$vc_coverage[d$vector_control==0] = 0
pop_size = 18.2 # in 100 thousands
npars = 4
serial_col = npars + 1
data_burnin = 5 # used 6 for timing plot
plot_years = 50  # used 5 for timing plot
last_col = serial_col + data_burnin + plot_years

tags = d[,1:npars]
data = d[,(serial_col + data_burnin + 1):last_col]
#medians = aggregate(d[,5:34], by=list(d$vector_control, d$vc_coverage, d$timing), FUN=median)

plot_effectiveness_over_time = function(tags, data, timing, plotcases=F, plotcumulative=F, plotcasesaverted=F) {
    #browser()
    medians = aggregate(data, by=list(tags$campaign_duration, tags$vc_coverage, tags$timing, tags$vector_control), FUN=median)
    names(medians)[1:npars]=c('dur','cov','day','vc')
    medians = medians[medians$day==timing | medians$vc==0,]
    eff_vals = 1 - sweep(as.matrix(medians[-1,(npars+1):(dim(medians)[2])]), 2, as.numeric(medians[1,(npars+1):(dim(medians)[2])]), '/')
    eff = cbind(medians[-1,1:npars], eff_vals)
    
    cum_medians = cbind(medians[,1:npars], t(apply(medians[,-c(1:npars)], 1, cumsum)))
    cum_eff_vals = 1 - sweep(as.matrix(cum_medians[-1,(npars+1):(dim(cum_medians)[2])]), 2, as.numeric(cum_medians[1,(npars+1):(dim(cum_medians)[2])]), '/')
    cum_eff = cbind(cum_medians[-1,1:npars], cum_eff_vals)
    reds=c('#FF0000','#DD0000','#AA0000')
    .lwd = 1
    duration = 1
    if (plotcases) {
       # browser()
        plotdata = rbind(medians[1,(npars+1):(npars+plot_years)], medians[medians$vc==1, (npars+1):(npars+plot_years)])/pop_size
        matplot(t(plotdata),
                lty=c(1,3:1), type='l', lwd=.lwd, col=c('black',reds),
                xlab='Year', ylab='Cases per 100,000 people', ylim=c(0,850), main='')
        legend('bottomright', legend=c('baseline','75% coverage','50% coverage','25% coverage'), lwd=.lwd, lty=c(1,1:3), bty='n', col=c('black',rev(reds)))
    } else if (plotcumulative) {
        matplot(t(cum_eff[cum_eff$dur==duration, (npars+1):(npars+plot_years)]),
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
        matplot(t(eff[eff$dur==duration,(npars+1):(npars+plot_years)]),
                lty=3:1, type='l', ylim=c(-0.25,1), lwd=.lwd, ylab='Effectiveness', xlab='Year', col=reds)
        abline(h=0,lty=2)
        legend('topright', legend=c('75% coverage','50% coverage','25% coverage'), lwd=.lwd, lty=1:3, bty='n', col=rev(reds))
    }
    return(medians)
}

#png('vc_effectiveness.png', width=1200, height=1920, res=180)
#png('vc_effectiveness-90_day_duration-jan1_start.png', width=3000, height=860, res=225)
#par(mfrow=c(1,3), mar=c(5.1,4.1,1,1), oma=c(0,0,2,0))
#
#par(las=1,bty='L')
#plot_effectiveness_over_time(tags,data,0, plotcases=T)
#mtext('Simulated impact of IRS (90-day campaign, January 1 start) on dengue cases per 100k, effectiveness, and cumulative effectiveness',side = 3,outer=T)
#plot_effectiveness_over_time(tags,data,0)
#plot_effectiveness_over_time(tags,data,0, plotcumulative=T)
#dev.off()

png('vc_4-panel-impact_50yr.png', width=2000, height=1720, res=225)
par(mfrow=c(2,2), mar=c(4.6,4.1,1.2,1), oma=c(0,0,2,0))

par(las=1,bty='L')
plot_effectiveness_over_time(tags,data,152, plotcases=T)
#grid()
#mtext('Simulated impact of IRS (90-day campaign, June 1 start)',side = 3,outer=T)
plot_effectiveness_over_time(tags,data,152)
#grid()
plot_effectiveness_over_time(tags,data,152, plotcasesaverted=T)
#grid()
plot_effectiveness_over_time(tags,data,152, plotcumulative=T)
#abline(v=1:10, h=1:20/20)
dev.off()

plot_years = 10
png('vc_4-panel-impact_10yr.png', width=2000, height=1720, res=225)
par(mfrow=c(2,2), mar=c(4.6,4.1,1.2,1), oma=c(0,0,2,0))

par(las=1,bty='L')
plot_effectiveness_over_time(tags,data[1:10],152, plotcases=T)
#grid()
#mtext('Simulated impact of IRS (90-day campaign, June 1 start)',side = 3,outer=T)
plot_effectiveness_over_time(tags,data[1:10],152)
#grid()
plot_effectiveness_over_time(tags,data[1:10],152, plotcasesaverted=T)
#grid()
plot_effectiveness_over_time(tags,data[1:10],152, plotcumulative=T)
#abline(v=1:10, h=1:20/20)
dev.off()
