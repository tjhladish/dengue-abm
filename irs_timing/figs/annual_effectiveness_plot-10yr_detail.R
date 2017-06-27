rm(list=ls())

require("RSQLite")
drv = dbDriver("SQLite")
db = dbConnect(drv, "./irs_stopping-effect_rerun.sqlite", flags=SQLITE_RO)

d <- dbGetQuery(db, 'select vector_control, timing, vc_coverage, campaign_duration, strat_years, M.*
                      from par P, met M, job J
                      where P.serial = M.serial
                    and strat_years = 50
                    and (vc_coverage = 0.75 or vector_control)
                    and P.serial = J.serial
                    and status = \'D\';')

pop_size = 18.2 # in 100 thousands
npars = 5
serial_col = npars + 1
data_burnin = 5 # used 6 for timing plot
plot_years = 10  # used 5 for timing plot
last_col = serial_col + data_burnin + plot_years

tags = d[,1:npars]
data = d[,(serial_col + data_burnin + 1):last_col]

plot_effectiveness_over_time = function(tags, data, timing) {
    medians = aggregate(data, by=list(tags$campaign_duration, tags$vc_coverage, tags$timing, tags$strat_years, tags$vector_control), FUN=median)
    names(medians)[1:npars]=c('dur','cov','day','strat','vc')
    medians = medians[medians$day==timing | medians$vc==0,]
    eff_vals = 1 - sweep(as.matrix(medians[-1,(npars+1):(dim(medians)[2])]), 2, as.numeric(medians[1,(npars+1):(dim(medians)[2])]), '/')
    eff = cbind(medians[-1,1:npars], eff_vals)
    reds=c('#FF0000','#DD0000','#AA0000')
    .lwd = 1:3
    matplot(t(eff[eff$dur==1 & eff$strat==50,(npars+1):(npars+plot_years)]),
            #lty=3:1, type='l', ylim=c(0,1.05), lwd=.lwd, ylab='Overall effectiveness', xlab='Year', col=reds)
            lty=3:1, type='l', lwd=.lwd, ylab='Overall effectiveness', xlab='Year', col=reds, ylim=c(0,1))
    #matlines(t(eff[eff$dur==1 & eff$strat==50, (npars+1):(npars+plot_years)]),
    #        lty=3:1, type='l', lwd=.lwd, main='', xlab='Year', col=reds)
    abline(h=0,lty=2)
    legend('topright', legend=c('75% coverage','50% coverage','25% coverage'), lwd=rev(.lwd), lty=1:3, bty='n', col=rev(reds))
    return(medians)
}

png('vc_effectiveness-main-refit-10year.png', width=1500, height=1290, res=200)
par(las=1,bty='L')
plot_effectiveness_over_time(tags,data,152)
#mtext('Simulated impact of IRS (90-day campaign, June 1) over first 10 years',side = 3, line=2)
dev.off()
