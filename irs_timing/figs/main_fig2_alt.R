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

plot_effectiveness_over_time = function(eff.dt, sero.dt, plot_years=21) {
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

    par(las=1)
    nf <- layout(matrix(c(1,2),ncol=1), widths=c(5,5,5), heights=c(3,2), TRUE)
#par(mar=c(1,4.2,2,1),oma=c(3.5,1,2,0))
    par(mar=c(1,4,0,0),oma=c(2.5,0,0,0))

#par(mfcol=c(2,1))

    .col=c('grey60', 'grey55', 'grey40', 'black') # increasing tone w/ increasing coverage
    .lwd = 3:6/2
    lty10=3; lty50=1 # ltys by 10 vs 50
    ylim_ = c(-4, 1.4)
    matplot(y=eff.ys, type='l',
            lwd=.lwd[-1], lty=c(rep(lty10,3), rep(lty50,3)),
            col=c(.col[-1], .col[-1]),
            ylab='Overall effectiveness', xlab='', ylim = ylim_, axes=F
           )
    lines(c(1,21), c(0,0), col=.col[1], lwd=.lwd[1])
    #abline(h=0,lty=3)
#    legend('bottomleft',
#           title=expression(bold("stuff")),
#           legend=c('75% coverage','50% coverage','25% coverage','75%, ended year 10','50%, ended year 10','25%, ended year 10'),
#           lwd=rev(.lwd), lty=1:3, col=c(rev(.col),rep('#999999',3)),
#           bty='n', cex=0.8)

    legend('bottomleft',
           inset=c(0.03,0.28),
           title=expression(bold("Intervention scenario")),
           title.adj=0.08,
           legend=c('Continued annual IRS campaigns', 'Campaigns ended after 10 years'),
           lwd=2.25, lty=c(lty50,lty10), col='black',
           bty='n', cex=0.8)

    legend('bottomleft',
           inset=c(0.0145,0.02),
           title=expression(bold("IRS household coverage")),
           title.adj=3.5,
           legend=c('75%','50%','25%', 'baseline (0%)'),
           lwd=rev(.lwd), lty=lty50, col=rev(.col),
           bty='n', cex=0.8)

    text(1, ylim_[2]-(0.025*(ylim_[2]-ylim_[1])), cex=2, font=2,labels='a')
    axis(1, at=0:4*5 + 1, labels=F)
    axis(2)

    matplot(y=sero.ys, type='l',
            lwd=c(.lwd[-1], .lwd),
            lty=c(rep(lty10,3), rep(lty50,4)),
            col=c(.col[-1], .col),
            ylab='Seroprevalence', xlab='', ylim = c(0,0.8), axes=F
           )

    ylim_ = par('usr')[3:4]
    text(1, ylim_[2]-(0.09*(ylim_[2]-ylim_[1])), cex=2, font=2,labels='b')
    axis(1, at=0:4*5 + 1, labels=F)
    axis(1, at=0:4*5 + 1, labels=0:4*5, lwd = 0, line = -0.3)
    axis(2)

#legend('bottomleft',
#       inset=c(0.025,0.05),
#       legend=c('Baseline (No IRS)','IRS, 75% coverage','IRS, 75% coverage ended in year 10'),
#       bty='n',
#       lwd=3,
#       col=.col[c(1,3,2)],
#       cex=0.8)

#avg_inf = t(mean1[rows,.range] + 2*mean2[rows,.range] + 3*mean3[rows,.range] + 4*mean4[rows,.range])
#matplot(avg_inf, type='l', col=.col, ylim=c(0,1.5), ylab='Mean number of past infections', xlab='',lwd=3,lty=1)
#legend('bottomright',inset=c(0.025,0.05),legend=c('Baseline (No IRS)','IRS, 75% coverage','IRS, 75% coverage ended in year 10'),bty='n',cex=1,lwd=3,col=.col[c(1,3,2)])
    mtext("Years since introducing IRS", side=1, outer=T, line=1, cex=1)


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
mag = 0.55*res/72
    png(args[3], width = 1000*mag, height = 1000*mag, units = "px", res=res)
    plot_effectiveness_over_time(effectiveness.dt, seroprevalence.dt)
dev.off()
