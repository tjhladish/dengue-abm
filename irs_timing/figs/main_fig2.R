rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)

# three cases
#  no intervention
#  75% coverage
#  75% coverage w/ stopping
effectiveness.dt  <- readRDS(args[1])

# should be 2 dimensions
# coverage (3 levels)
# stopping vs not (2 levels)
seroprevalence.dt <- readRDS(args[2])

# two facets:
## effectiveness
## seroprevalence

plot_years = 20  # used 5 for timing plot
# last_col = serial_col + data_burnin + plot_years
# 
# tags = d[,1:npars]
# data = d[,(serial_col + data_burnin + 1):last_col]
# 
# # data0 = d[,60:84]    ## UGLY HACK
# # data1 = d[,113:137]  ## UGLY HACK
# # data2 = d[,166:190]  ## UGLY HACK
# 
# data0 = d[,62:86]    ## UGLY HACK
# data1 = d[,117:141]  ## UGLY HACK
# data2 = d[,172:196]  ## UGLY HACK
# data3 = d[,227:251]  ## UGLY HACK
# data4 = d[,282:306]  ## UGLY HACK
# 
# mean0 = aggregate(data0, by=list(tags$campaign_duration, tags$vc_coverage, tags$timing, tags$strat_years, tags$vector_control), FUN=mean)
# mean1 = aggregate(data1, by=list(tags$campaign_duration, tags$vc_coverage, tags$timing, tags$strat_years, tags$vector_control), FUN=mean)
# mean2 = aggregate(data2, by=list(tags$campaign_duration, tags$vc_coverage, tags$timing, tags$strat_years, tags$vector_control), FUN=mean)
# mean3 = aggregate(data3, by=list(tags$campaign_duration, tags$vc_coverage, tags$timing, tags$strat_years, tags$vector_control), FUN=mean)
# mean4 = aggregate(data4, by=list(tags$campaign_duration, tags$vc_coverage, tags$timing, tags$strat_years, tags$vector_control), FUN=mean)

.col = c('#cce3ff','#99c7ff','#4a91cc','#2f5b7f','#104060')
stackimmunity = function(r,.label='') {
    .range = 11:30
    mat = matrix(c(mean0[r,.range],mean1[r,.range],mean2[r,.range],mean3[r,.range],mean4[r,.range]), ncol=5, nrow=length(.range))
    stackpoly(mat, stack=T, ylim=c(0,1), col=.col, axis4=F, xaxlab = F, axis=F)
    axis(2,tick = T,labels=c('0.0','0.2','0.4','0.6','0.8','1.0'),at=0:5/5, las=1)
    axis(4,labels = F)
    text(1.2, 0.07, pos=4, labels=.label, font=2)
}

#legend('topright',legend = c('all houses treated in 1 day', 'houses treated across 90 days', 'houses treated continuously'), lty=c(2,1,3), lwd=lwd_, bty='n')

png('susceptible_and_mean_immunity.png',height=1200,width=1000,res=180)
#mean0 = aggregate(data0, by=list(tags$campaign_duration, tags$vc_coverage, tags$timing, tags$strat_years, tags$vector_control), FUN=mean)
rows = c(1,7,4)
.range = 11:30
par(mfrow=c(2,1),mar=c(1,4.2,2,1),oma=c(3.5,1,0,0))
matplot(t(mean0[rows,.range]), type='l', col='black', ylim=c(0,0.6), ylab='Fully susceptible population', xlab='')
legend('bottomright',inset=c(0.025,0.05),legend=c('Baseline (No IRS)','IRS, 75% coverage','IRS, 75% coverage ended in year 10'), lty=1:3,bty='n',cex=0.75)
avg_inf = t(mean1[rows,.range] + 2*mean2[rows,.range] + 3*mean3[rows,.range] + 4*mean4[rows,.range])
matplot(avg_inf, type='l', col='black', ylim=c(0,1.5), ylab='Mean number of past infections', xlab='')
mtext("Intervention year", side=1, outer=T, line=1.5)
dev.off()


plot_effectiveness_over_time = function(tags, data, timing) {
    medians = aggregate(data, by=list(tags$campaign_duration, tags$vc_coverage, tags$timing, tags$strat_years, tags$vector_control), FUN=median)
    names(medians)[1:npars]=c('dur','cov','day','strat','vc')
    medians = medians[medians$day==timing | medians$vc==0,]
    eff_vals = 1 - sweep(as.matrix(medians[-1,(npars+1):(dim(medians)[2])]), 2, as.numeric(medians[1,(npars+1):(dim(medians)[2])]), '/')
    eff = cbind(medians[-1,1:npars], eff_vals)
    reds=c('#FF0000','#DD0000','#AA0000')
    .lwd = 1:3
    ylim_ = c(-4.4, 1.4)
    matplot(t(eff[eff$dur==1 & eff$strat==10,(npars+1):(npars+plot_years)]),
            #lty=3:1, type='l', ylim=c(0,1.05), lwd=.lwd, ylab='Overall effectiveness', xlab='Year', col=reds)
            lty=3:1, type='l', lwd=.lwd, ylab='Overall effectiveness', xlab='', col='#999999', ylim=ylim_, axes=F)
    matlines(t(eff[eff$dur==1 & eff$strat==50, (npars+1):(npars+plot_years)]),
             lty=3:1, type='l', ylab='Cumulative cases averted per 100 people', lwd=.lwd, col=reds)
    #browser()
    #ylim_ = par('usr')[3:4]
    text(1, ylim_[2]-(0.025*(ylim_[2]-ylim_[1])), cex=2, font=2,labels='a')
    box()
    axis(1, labels=F)
    axis(1, lwd = 0, line = -0.3)
    axis(2)

    abline(h=0,lty=2)
    legend('bottomleft', legend=c('75% coverage','50% coverage','25% coverage','75%, ended year 10','50%, ended year 10','25%, ended year 10'), lwd=rev(.lwd), lty=1:3, bty='n', col=c(rev(reds),rep('#999999',3)),cex=0.8)
    return(medians)
}

png('vc_effectiveness-w-immunity-20year-main.png', width=1500, height=1500, res=200)
#par(las=1,bty='L')
par(las=1)
nf <- layout(matrix(c(1,2),ncol=1), widths=c(5,5,5), heights=c(3,2), TRUE)
#par(mar=c(1,4.2,2,1),oma=c(3.5,1,2,0))
par(mar=c(1,4,1,0),oma=c(3,0,0,0))
plot_effectiveness_over_time(tags,data,152)
#mtext('Simulated impact of IRS (90-day campaign, June 1 start) over 10-20 years',side = 3, line=2,cex=0.8)

rows = c(1,4,7)
.range = 11:30
.col = c('black','#999999','#AA0000')
matplot(1-t(mean0[rows,.range]), type='l', col=.col, ylim=c(0,0.8), ylab='Seroprevalence', xlab='',lwd=3,lty=1, axes=F)
ylim_ = par('usr')[3:4]
text(1, ylim_[2]-(0.09*(ylim_[2]-ylim_[1])), cex=2, font=2,labels='b')
box()
axis(1, labels=F)
axis(1, lwd = 0, line = -0.3)
axis(2)

legend('bottomleft',inset=c(0.025,0.05),legend=c('Baseline (No IRS)','IRS, 75% coverage','IRS, 75% coverage ended in year 10'),bty='n',lwd=3,col=.col[c(1,3,2)],cex=0.8)
avg_inf = t(mean1[rows,.range] + 2*mean2[rows,.range] + 3*mean3[rows,.range] + 4*mean4[rows,.range])
#matplot(avg_inf, type='l', col=.col, ylim=c(0,1.5), ylab='Mean number of past infections', xlab='',lwd=3,lty=1)
#legend('bottomright',inset=c(0.025,0.05),legend=c('Baseline (No IRS)','IRS, 75% coverage','IRS, 75% coverage ended in year 10'),bty='n',cex=1,lwd=3,col=.col[c(1,3,2)])
mtext("Intervention year", side=1, outer=T, line=1, cex=1)
dev.off()
