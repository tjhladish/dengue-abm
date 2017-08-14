rm(list=ls())

require(data.table)

args <- commandArgs(trailingOnly = T)
# args <- paste0("~/Dropbox/who/fig1_data/", c("baseline.cases.rds", "interventions.cases.rds", "averted.cases.rds", "stopping-eff.rds", "vc_4-panel-impact_10yr.png"))

plot_years <- as.integer(gsub(".+_(\\d+)yr.png","\\1", args[5]))

baseline.cases.dt <- readRDS(args[1])[year < plot_years,.(
  cases.md = median(cases), cum.cases.md = median(cum.cases)
), by=year]
interventions.cases.dt <- readRDS(args[2])[year < plot_years,.(
  cases.md = median(cases), cum.cases.md = median(cum.cases)
), by=.(coverage, year)]
averted.cases.dt <- readRDS(args[3])[year < plot_years]
eff.dt <- readRDS(args[4])[end_year == 50 & year < plot_years]

pop_size = 18.2 # in 100 thousands

bcol <- 'black'
cols=c('grey55', 'grey40', 'black') # increasing tone w/ increasing coverage

plot_effectiveness_over_time = function(
  bcases, icases, acases, eff,
  plotcases=F, plotcumulative=F, plotcasesaverted=F
) {
    reds=cols #c('#FF0000','#DD0000','#AA0000')
    .lwd = 1
    duration = 1
    if (plotcases) {
       # browser()
        plotdata = cbind(bcases$cases.md, matrix(icases[,cases.md,keyby=.(coverage,year)]$cases.md, byrow = F, ncol=3))/pop_size
        matplot(plotdata,
                lty=c(3,rep(1,3)), type='l', lwd=.lwd, col=c('black',reds),
                xlab='Year', ylab='Cases per 100,000 people', ylim=c(0,850), main='')
        legend('bottomright', legend=c('baseline','75% coverage','50% coverage','25% coverage'), lwd=.lwd, lty=c(3,rep(1,3)), bty='n', col=c('black',rev(reds)))
    } else if (plotcumulative) {
        plotdata = matrix(acases[,cum.eff.md,keyby=.(coverage,year)]$cum.eff.md, byrow = F, ncol=3)
        matplot(plotdata,
                lty=1, type='l', ylim=c(0,1), ylab='Cumulative effectiveness', lwd=.lwd, main='', xlab='Year', col=reds)
        legend('topright', legend=c('75% coverage','50% coverage','25% coverage'), lwd=.lwd, lty=1, bty='n', col=rev(reds))
    } else if (plotcasesaverted) {
        plotdata = matrix(acases[,cum.averted.md,keyby=.(coverage,year)]$cum.averted.md, byrow = F, ncol=3)/(pop_size*1000)
        matplot(plotdata,
                lty=1, type='l', ylab='Cumulative cases averted per 100 people', lwd=.lwd, main='', xlab='Year', col=reds)
        legend('bottomright', legend=c('75% coverage','50% coverage','25% coverage'), lwd=.lwd, lty=1, bty='n', col=rev(reds))
    } else {
        plotdata = matrix(eff[,q.med,keyby=.(coverage,year)]$q.med, byrow = F, ncol=3)
        matplot(plotdata,
                lty=1, type='l', ylim=c(-0.25,1), lwd=.lwd, ylab='Effectiveness', xlab='Year', col=reds)
        abline(h=0,lty=2)
        legend('topright', legend=c('75% coverage','50% coverage','25% coverage'), lwd=.lwd, lty=1, bty='n', col=rev(reds))
    }
}

png(args[5], width=2000, height=1720, res=225)
par(mfrow=c(2,2), mar=c(4.6,4.1,1.2,1), oma=c(0,0,2,0))

par(las=1,bty='L')
plot_effectiveness_over_time(baseline.cases.dt, interventions.cases.dt, averted.cases.dt, eff.dt, plotcases=T)
#grid()
#mtext('Simulated impact of IRS (90-day campaign, June 1 start)',side = 3,outer=T)
plot_effectiveness_over_time(baseline.cases.dt, interventions.cases.dt, averted.cases.dt, eff.dt)
#grid()
plot_effectiveness_over_time(baseline.cases.dt, interventions.cases.dt, averted.cases.dt, eff.dt, plotcasesaverted=T)
#grid()
plot_effectiveness_over_time(baseline.cases.dt, interventions.cases.dt, averted.cases.dt, eff.dt, plotcumulative = T)
#abline(v=1:10, h=1:20/20)
dev.off()
