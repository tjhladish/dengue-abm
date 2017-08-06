rm(list=ls())

require(data.table)
require(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
# args <- c(paste0("~/Dropbox/who/fig1_data/stopping-",c("eff","sero"),".rds"),"~/Dropbox/who/fig1_data/fig2_alt.png")

# should be 2 dimensions
# coverage (3 levels)
# stopping vs not (2 levels)
effectiveness.dt  <- readRDS(args[1])

# three cases
#  no intervention
#  75% coverage
#  75% coverage w/ stopping
seroprevalence.dt <- readRDS(args[2])

# two facets:
## effectiveness
## seroprevalence

eff.slice.dt <- effectiveness.dt[, .(value=q.med,measure="Annual effectiveness"), keyby=.(end_year, coverage, Year=year)]
sero.slice.dt <- seroprevalence.dt[, .(value=q.med,measure="Seroprevalence"), keyby=.(end_year, coverage, Year=year)]
seroprevalence.dt[end_year == 0, end_year := 50L]

plot_effectiveness_over_time = function(eff.dt, sero.dt, plot_years=20) {
  slice.eff.dt  <- eff.dt[year < plot_years, q.med, keyby=.(end_year, coverage, year)]
  slice.sero.dt <- sero.dt[year < plot_years, q.med, keyby=.(end_year, coverage, year)]
  # keyby sets order ascending by end_year, coverage, year
  # increasing 'end_year' -> all the 10s before the 50s
  # next, increasing 'coverage' -> within 10/50, 25 -> 50 -> 75
  # last, increasing year w/in each scenario
  # this is going to lead to drawing 10 first, then 50

  # build matrix column-by-column from slice.dt
  eff.ys <- matrix(slice.eff.dt$q.med, byrow=F, nrow=plot_years)
  sero.ys <- matrix(slice.sero.dt$q.med, byrow=F, nrow=plot_years)
  
  curpar <- par()
  par(mfcol=c(2,1))
  
  reds=c('lightgrey','grey','black') # increasing tone w/ increasing coverage
  .lwd = 1 # only one line weight
  lty10=3; lty50=1 # ltys by 10 vs 50
  ylim_ = c(-4.4, 1.4)
  matplot(y=eff.ys, type='l',
    lwd=.lwd, lty=c(rep(lty10,3), rep(lty50,3)),
    col=c(reds, reds),
    ylab='Overall effectiveness', xlab='', ylim = ylim_ #, axes=F
  )
  
  matplot(y=sero.ys, type='l',
          lwd=.lwd, lty=c(rep(lty10,3), rep(lty50,4)),
          col=c(reds,'grey80', reds),
          ylab='Seroprevalence', xlab='', ylim = c(0,0.8) #, axes=F
  )

  ## tjh: TODO, if you want legend-y stuff sorted
  # text(1, ylim_[2]-(0.025*(ylim_[2]-ylim_[1])), cex=2, font=2,labels='a')
  # box()
  # axis(1, labels=F)
  # axis(1, lwd = 0, line = -0.3)
  # axis(2)
  # 
  # abline(h=0,lty=2)
  # legend('bottomleft', legend=c('75% coverage','50% coverage','25% coverage','75%, ended year 10','50%, ended year 10','25%, ended year 10'), lwd=rev(.lwd), lty=1:3, bty='n', col=c(rev(reds),rep('#999999',3)),cex=0.8)
  par(curpar)
}

res = 300
mag = 0.85*res/72
png(args[3], width = 1000*mag, height = 1000*mag, units = "px", res=res)
plot_effectiveness_over_time(effectiveness.dt, seroprevalence.dt)
dev.off()
