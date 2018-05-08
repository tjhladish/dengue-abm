rm(list=ls())

args <- c(
  "~/Dropbox/who/fig1_data/stat.eff10.rds",
  "~/Dropbox/who/fig1_data/extra-stat.eff10.rds",
  "~/Dropbox/who/fig1_data/R0.rds",
  "~/Dropbox/who/fig1_data/mos.rds",
  "~Dropbox/who/fig1_data/irs_timing-realigned-SI.png"
)

args <- commandArgs(trailingOnly = T)

require(data.table)

stat.eff.dt <- rbind(
  readRDS(args[1]), readRDS(args[2])
)

R0.dt <- readRDS(args[3])
mos.dt <- readRDS(args[4])

month_starts = c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)
month_labels = c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec')

# wrap days
stat.eff.dt[
  duration != 365 & durability == 90 & efficacy == 0.8,
  moddoy := doy + floor((duration+durability)/2)
]
stat.eff.dt[moddoy > 365, moddoy := moddoy - 365 ]


backcols <- c(mos="blue",r0=3)
# `Mos. pop.`="blue",
# `R0`=3,

png(tail(args,1), width=2000, height=2200, res=300)
par(mfrow=c(3,1), mar=c(1,2.1,1,1), oma=c(3,3,0,0),las=1)
#dur_ = 1
end = F
for (cov_ in c(75,50,25)) {
    test = stat.eff.dt[coverage==cov_ & durability == 90 & efficacy == 0.8, .(med.eff10, smooth, doy), keyby=.(duration, moddoy)]
    minval = test[, min(med.eff10)]
    maxval = test[, max(med.eff10)]
    .range = maxval-minval
    pad    = 0.05
    plot(0:365, ylim=c(minval-pad*.range,maxval+pad*.range), xlab='', ylab='', type='n', xaxt='n')
    axis(1,at=c(month_starts,365), labels = F) #,labels = month_labels)
    axis(1,at=month_starts+15, tick = F, labels = month_labels, padj = -1.5)
    
    with(R0.dt[layer == "foreground" & duration != 365], lines(doy, (.range*value)+minval, col=backcols["r0"], lwd=2)) # R0
    with(mos.dt[layer == "foreground"], lines(doy, (.range*value)+minval, col=backcols["mos"], lwd=1.5)) # mosquitoes

    # shift ts based on half campaign duration + half durability

    with(test[duration == 1],{
      doyorder <- order(doy)
      lines(doy[doyorder], med.eff10[doyorder], type='l', lty=2, lwd=2, col="grey")
      #lines(moddoy, smooth, type='l', lty=1, lwd=2)
    })
    
    with(test[duration == 90],{
      doyorder <- order(doy)
      lines(doy[doyorder], med.eff10[doyorder], type='l', lty=1, lwd=2, col="grey")
      #lines(moddoy, smooth, type='l', lty=2, lwd=2)
    })
    
        
    with(test[duration == 1],{
      doyorder <- order(doy)
      lines(moddoy, med.eff10, type='l', lty=2, lwd=2)
      #lines(moddoy, smooth, type='l', lty=1, lwd=2)
    })
    
    with(test[duration == 90],{
      doyorder <- order(doy)
      lines(moddoy, med.eff10, type='l', lty=1, lwd=2)
      #lines(moddoy, smooth, type='l', lty=2, lwd=2)
    })

    with(test[duration == 365], abline(h=med.eff10, lty=3, lwd=2))
    
    #legend('bottomright',inset = c(0.08,0.16), legend = paste0(cov_*100, '% household coverage'), bty='n')
    text(160, minval, pos=4, labels = paste0(cov_, '% household coverage'), font=2)
    legend('topleft',title='Rollout period', legend = c('1 day', '90 days', '365 days'), lty=c(2,1,3), bty='n', lwd=2)
}
par(las=0)
mtext("Timing of peak coverage (Month)", side=1, outer=T, line=1)
mtext("10 year effectiveness (prevented cases)", side=2, outer=T, line=1)
#mtext("Effect of peak insecticide coverage on IRS effectiveness", side=3, outer=T)
dev.off()
