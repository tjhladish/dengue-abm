rm(list=ls())

require(data.table)

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

    par(mfrow=c(3,1))
    #nf <- layout(matrix(c(1,2),ncol=1), widths=c(5,5,5), heights=c(3,2), TRUE)
#par(mar=c(1,4.2,2,1),oma=c(3.5,1,2,0))
    par(las=1,mar=c(1,5,0,0),oma=c(2.5,0,0,0))

#par(mfcol=c(2,1))

    #col50='black'
    col0 ='black'
    col10='#aa0000'
    col50=c('grey55', 'grey40', 'black') # increasing tone w/ increasing coverage
    .lwd = 3
    #.lwd = c(2,1.5,2.25,3)
    #.lwd = 3:6/2
    lty0 =3
    lty10=1
    lty50=1
    ylim_ = c(0, 1.1)

    #### PANEL A -- Long-term annual effectiveness
    matplot(y=cbind(rep(0,21),eff.ys[,4:6]), type='l',
            lwd=.lwd, lty=c(lty0,rep(lty50,3)),
            col=c(col0, rep(col50,3)),
            #col=c(.col[-1], .col[-1]),
            ylab='Overall effectiveness', xlab='', ylim = ylim_, axes=F,
            cex.lab=1.5
           )
    #lines(c(1,21), c(0,0), col=col50, lwd=.lwd[1])
    text(1, ylim_[2]-(0.025*(ylim_[2]-ylim_[1])), cex=2, font=2,labels='a')
    axis(1, at=0:4*5 + 1, labels=F)
    axis(2)

    legend('topright',
           inset=c(0.05,0.05),
           legend=c('75% coverage','50% coverage','25% coverage', 'baseline (0% coverage)'),
           lwd=rev(.lwd), lty=c(rep(lty50,3),lty0), col=rev(c(col0, rep(col50,3))),
           seg.len=2.5,
           #lwd=rev(.lwd), lty=lty50, col=rev(.col),
           bty='n', cex=1.25)


    #### PANEL B -- Long-term annual effectiveness w/ stopping
    ylim_ = c(-4, 1.5)
    matplot(y=cbind(rep(0,21),eff.ys[,c(3,6)]), type='l',
            lwd=.lwd, lty=c(lty0,lty10,lty50),
            col=c(col0, col10, col50[3]),
            #col=c(.col[-1], .col[-1]),
            ylab='Overall effectiveness', xlab='', ylim = ylim_, axes=F,
            cex.lab=1.5
           )

    legend('bottomleft',
           inset=c(0.05,0.05),
           legend=c('75% coverage','75% coverage ending in year 10','baseline (0% coverage)'),
           lwd=rev(.lwd), lty=c(lty50,lty10,lty0), col=c(col50[3],col10,col0),
           seg.len=2.5,
           #lwd=rev(.lwd), lty=lty50, col=rev(.col),
           bty='n', cex=1.25)

    text(1, ylim_[2]-(0.025*(ylim_[2]-ylim_[1])), cex=2, font=2,labels='b')
    axis(1, at=0:4*5 + 1, labels=F)
    axis(2)

    #### PANEL C -- Long-term annual seroprevalence
    matplot(y=sero.ys[,c(3,4,7)], type='l',
            lwd=.lwd,
            lty=c(lty10,lty0,lty50),
            col=c(col10, col0, col50[3]),
            ylab='Seroprevalence', xlab='', ylim = c(0.3,0.8), axes=F,
            cex.lab=1.5
           )

    legend('bottomleft',
           inset=c(0.05,0.05),
           legend=c('75% coverage','75% coverage ending in year 10','baseline (0% coverage)'),
           lwd=rev(.lwd), lty=c(lty50,lty10,lty0), col=c(col50[3],col10,col0),
           seg.len=2.5,
           #lwd=rev(.lwd), lty=lty50, col=rev(.col),
           bty='n', cex=1.25)

    ylim_ = par('usr')[3:4]
    text(1, ylim_[2]-(0.09*(ylim_[2]-ylim_[1])), cex=2, font=2,labels='c')
    axis(1, at=0:4*5 + 1, labels=F)
    axis(1, at=0:4*5 + 1, labels=0:4*5, lwd = 0, line = -0.3)
    axis(2)

    mtext("Years since introducing IRS", side=1, outer=T, line=1, cex=1)

    par(curpar)
}

res = 300
mag = 0.55*res/72
    png(args[3], width = 1000*mag, height = 1200*mag, units = "px", res=res)
    plot_effectiveness_over_time(effectiveness.dt, seroprevalence.dt)
dev.off()
