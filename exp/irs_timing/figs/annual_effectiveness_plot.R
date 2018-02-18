rm(list=ls())

require(data.table)

args <- commandArgs(trailingOnly = T)
# args <- c(paste0("~/Dropbox/who/fig1_data/stopping-", c("baseline.cases.rds", "interventions.cases.rds", "averted.rds", "eff.rds")), "~/Dropbox/who/fig1_data/vc_4-panel-impact_10yr.png")
# args <- c("~/Dropbox/who/fig1_data/foi-baseline.cases.rds", "~/Dropbox/who/fig1_data/foi-interventions.cases.rds", "~/Dropbox/who/fig1_data/foi-averted.rds", "~/Dropbox/who/fig1_data/foi-eff.rds", "~/Dropbox/who/fig1_data/foi_4panel_20yr.png")

plot_years <- as.integer(gsub(".+_(\\d+)yr.png","\\1", args[5]))

base.dt <- readRDS(args[1])
bkeys <- grep("(case|particle)", names(base.dt), invert=T, value=T)

baseline.cases.dt <- base.dt[year < plot_years,.(
  cases.md = median(cases), cum.cases.md = median(cum.cases)
), by=bkeys ]

i.dt <- readRDS(args[2])
ikeys <- grep("(case|particle)", names(i.dt), invert=T, value=T)
plotkeys <- grep("year",ikeys,invert = T, value = T)

interventions.cases.dt <- i.dt[end_year == 50 & year < plot_years,.(
  cases.md = median(cases), cum.cases.md = median(cum.cases)
), by=ikeys]

averted.cases.dt <- readRDS(args[3])[end_year == 50 & year < plot_years]
eff.dt <- readRDS(args[4])[end_year == 50 & year < plot_years]

pop_size = 18.2 # in 100 thousands

bcol <- 'black'
cols=c('grey55', 'grey40', 'black') # increasing tone w/ increasing coverage

plot_effectiveness_over_time = function(
  bcases, icases, acases, eff,
  plotcases=F, plotcumulative=F, plotcasesaverted=F
) {
    reds=cols #c('#FF0000','#DD0000','#AA0000')
    .lwd = 2
    duration = 1
    if (plotcases) {
       # browser()
        plotdata = cbind(bcases$cases.md, matrix(icases[,cases.md,keyby=.(coverage,year)]$cases.md, byrow = F, ncol=3))/pop_size
        matplot(plotdata,
                lty=c(3,rep(1,3)), type='l', lwd=.lwd, col=c('black',reds),
                xlab='Year', ylab='Cases per 100,000 people', ylim=c(0,1220), main='')
        legend('bottomright', legend=c('baseline','75% coverage','50% coverage','25% coverage'), seg.len=1.5, lwd=.lwd, lty=c(3,rep(1,3)), bty='n', col=c('black',rev(reds)))
    } else if (plotcumulative) {
        plotdata = matrix(acases[,cum.eff.md,keyby=.(coverage,year)]$cum.eff.md, byrow = F, ncol=3)
        matplot(plotdata,
                lty=1, type='l', ylim=c(0,1), ylab='Cumulative effectiveness', lwd=.lwd, main='', xlab='Year', col=reds)
        legend('topright', legend=c('75% coverage','50% coverage','25% coverage'), seg.len=1.5, lwd=.lwd, lty=1, bty='n', col=rev(reds))
    } else if (plotcasesaverted) {
        plotdata = matrix(acases[,cum.averted.md,keyby=.(coverage,year)]$cum.averted.md, byrow = F, ncol=3)/(pop_size*1000)
        matplot(plotdata,
                lty=1, type='l', ylab='Cumulative cases averted per 100 people', lwd=.lwd, main='', xlab='Year', col=reds)
        legend('bottomright', legend=c('75% coverage','50% coverage','25% coverage'), seg.len=1.5, lwd=.lwd, lty=1, bty='n', col=rev(reds))
    } else {
        plotdata = matrix(eff[,q.med,keyby=.(coverage,year)]$q.med, byrow = F, ncol=3)
        matplot(plotdata,
                lty=1, type='l', ylim=c(-0.25,1), lwd=.lwd, ylab='Effectiveness', xlab='Year', col=reds)
        abline(h=0,lty=2)
        legend('topright', legend=c('75% coverage','50% coverage','25% coverage'), seg.len=1.5, lwd=.lwd, lty=1, bty='n', col=rev(reds))
    }
}

png(args[5], width=2000, height=1720, res=225)
par(mfrow=c(2,2), mar=c(4.6,4.1,1.2,1), oma=c(0,0,2,0), lend=2)

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
