rm(list=ls())
library(scales)
require(plotrix)
require(data.table)
julian_days = 1:365
# days on which to plot weekly, monthly data
week_days = c(1, seq(3.5,365,7), 365)
# TODO - make sure that ticks and day-of-month stuff is coherent
month_days = c(31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365) - 15
month_days2 = c(1,month_days,365)
dom = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366)
months = format(ISOdatetime(2000,1:12,1,0,0,0),"%b")
months = c(months, '') # "" for Jan at end of year
orange = '#ff7700'

########################################
##### Panel 1: Temperature --> EIP
########################################
e = read.table('seasonal_eip.bad', header=T)
e_avg = read.table('seasonal_EIP_24hr.out', header=F)
max_eip = max(e_avg[,1])
rescaled_e_avg = e_avg[,1]/max_eip
t = read.table('merida_temps.out', header=F)
max_y1 = max(t[,1])
rescaled_tmin = e[,3]/max_y1
rescaled_tmax = e[,4]/max_y1
rescaled_t    = t[,1]/max_y1

########################################
##### Panel 2: Rain --> Mosquito pop
########################################
all_data = read.csv("NOAA_yucatan_daily.csv")

# work with Merida data only
d = all_data[all_data$NAME=='AEROP.INTERNACIONAL',]
d$DATE = as.Date(d$DATE)

precip=data.frame(DATE=seq(as.Date("1979/1/1"), as.Date("2014/1/1"), "days"))
precip$PRCP = d$PRCP[match(precip$DATE, d$DATE)]

precip$bool0 = precip$PRCP > 0
precip$month = strftime(precip$DATE, '%m')
precip$day = strftime(precip$DATE, '%m-%d')
seasonal_rain = aggregate(precip$bool0, by=list(precip$day), mean, na.rm=T)

names(seasonal_rain) = c('day', 'precip')
seasonal_rain = seasonal_rain[-which(seasonal_rain$day == '02-29'),]
mod3 = smooth.spline(rep(seasonal_rain$precip, 3),spar=0.6)
seasonal_rain$smooth = mod3$y[1:365 + 365]

write.table(seasonal_rain[, c('smooth', 'precip', 'day')], 'merida_preciptiation.out', row.names=F, col.names=F, quote=F)

seasonal_rain$mosquitoes = seasonal_rain$smooth/max(seasonal_rain$smooth)
# lag by one week to account for egg-to-emergence time, based on Table 2 in Rueda et al, Temperature-dependent
# development and survival rates of Culex quinquefasciatus and Aedes aegypti (Diptera: Culicidae), J Med Ent, 1990, 27:5
seasonal_rain$mosquitoes = seasonal_rain$mosquitoes[(0:364-7)%%365 + 1]

write.table(seasonal_rain[, 'mosquitoes'], 'mosquito_seasonality.out', row.names=F, col.names=F, quote=F)

max_rain = max(seasonal_rain$precip)

pdf('rain_fit.pdf', width=10, height=6)
plot(seasonal_rain$precip, type='l', lwd=0.5)
lines(seasonal_rain$smooth,col=2, lwd=2)
lines(1:365+8, seasonal_rain$smooth,col=3, lwd=2)
dev.off()

########################################
##### Panel 3: R0
########################################
library("RSQLite")
drv = dbDriver("SQLite")
db = dbConnect(drv, "rzero-refit.sqlite", flags=SQLITE_RO)
all_r0 = dbGetQuery(db, 'select realization, m.* from met m, par p where m.serial=p.serial')
r0_100 = aggregate(all_r0, by=list(all_r0$realization), FUN=mean)
r0_100 = r0_100[,-c(1:3)]
r0 = sapply(r0_100, mean)
r0_q = sapply(r0_100, quantile, probs=c(0.025, 0.975))
r0_q025 = t(r0_q)[,1]
r0_q975 = t(r0_q)[,2]
smooth_tmp = smooth.spline(rep(r0, 3),spar=0.5)
r0_smooth = smooth_tmp$y[1:365 + 365]
#plot(r0_smooth, type='l')
#abline(h=1)
r0_gt1 = which(r0_smooth>=1)
season_start = r0_gt1[1]
season_end   = r0_gt1[length(r0_gt1)]

########################################
##### Panel 4: Summaries + Epi data
########################################
#reported = read.table('aggregated_cases.out', header=T)
#epi = aggregate(reported$cases, by=list(reported$week), mean)$x

obs = read.table('weekly_confirmed_cases_1995-2015_wide.csv', header=F, sep=' ', col.names=1:53, fill=T)
epi = colMeans(obs[,1:52], na.rm=T)

#week_days = 0:51*7+3.5
epi_jan1 = epi[1]*(4.5/8) + epi[52]*(3.5/8) # weight according to inverse distance from jan1/dec31
epi = c(epi_jan1, epi, epi_jan1)
epi = epi/max(epi)

# simulation data
# tag day intros local infec mild severe eip nmos ef_mild ef_severe
# output file digested with:
# tail -n +2 posterior_daily_yucatan_1995-2011.out | awk '{print ($2+100)%365, $6/$10, $7/$11}' > posterior_daily_yucatan_1995-2011.digest
# ./aggregate_years.py

# newer version, reflecting changes to log format:
# grep -h day abc_yucatan-sm3.err-* | awk '{if ($4 >= 42240 && $4 < 48445) print $1, $4, ($4+100)%365, $8*$11, $9*$12}' > ../posterior_daily_1995-2011-sm3.alt
# sort --parallel 4 posterior_daily_1995-2011-sm3.alt | uniq | awk '{print $3, $4, $5}' > posterior_daily_1995-2011-sm3.digest
#sim_epi = read.table('sim_days', header=F)

load(file='sim.dt.mean.rda')
sim.overall.mean = aggregate(sim.dt.mean[,.(inc)], by=list(day=sim.dt.mean$day), FUN=mean)
sim_epi = sim.overall.mean$inc


plot_seasonality = function(supplement=T) {
    filename = 'seasonality_detail.png'
    hres = 1750
    wres = 1500
    res  =  240
    if (!supplement) {
        filename = 'seasonality_summary.png'
        hres = 1000
        res  = 200
    }

    png(filename, width=wres, height=hres, res=res);
    par(las=1)
    par(mfrow=c(2,1))
    par(mar=c(2.0, 4.1, 0, 4.1))
    par(oma=c(2.5,0.5,0.5,0.5))

    if (supplement) {
        par(mfrow=c(4,1))
        # temp / EIP panel
        plot(rescaled_t, col=1, lwd=1.5, xlab='', type='l',ylim=c(0,1.5), axes=F, ylab=expression(paste('Temperature (',degree,'C)')))
        usr.lim = par('usr')
        box(bty='c')
        lines(rescaled_tmin, col=1, lwd=0.5)
        lines(rescaled_tmax, col=1, lwd=0.5)
        lines(rescaled_e_avg, col=3, lwd=1.5)
        axis(1, at=dom, labels=rep('',13))
        axis(2, at=0:3*15/max_y1, labels=0:3*15)
        axis(4, at=usr.lim[3:4], labels=F, col=3, lwd.tick=0)
        axis(4, at=0:2*10/max_eip, labels=0:2*10, col=3, col.axis=3, lwd=0, lwd.tick=1)
        mtext("Incubation period (days)", 4, line=3, cex=0.75, las=3, col=3)
        legend('topright', legend=c('Min/max temp', 'Mean temp', 'EIP'), col=c(1,1,3), lwd=c(0.5,1.5,1.5), bty='n', cex=0.8)

        # precip / mos panel
        max_y1 = max(seasonal_rain$smooth)
        max_y2 = max(seasonal_rain$mosquitoes)
        plot(seasonal_rain$precip/max_y1, type='l', lwd=0.5, xlab='', axes=F, ylab='Pr{precipitation}')
        print(par('usr'))
        box(bty='c')
        lines(seasonal_rain$smooth/max_y1, col=1, lwd=1.5);
        lines(1:365, seasonal_rain$mosquitoes/max_y2,col=4, lwd=1.5)
        legend('topright', legend=c('Pr{precipitation}', 'Smoothed precip', 'Mosquito population'), col=c(1,1,4), lwd=c(0.5,1.5,1.5), lty=c(1,1,1), bty='n', cex=0.8)
        axis(1, at=dom, labels=rep('',13))
        axis(2, at=((0:3)*0.2)/max_y1, labels=(0:3*0.2))
        axis(4, at=usr.lim[3:4], labels=F, col=4, lwd.tick=0)
        axis(4, col=4, col.axis=4, lwd=0, lwd.tick=1)
        mtext("Relative population size", 4, line=3, cex=0.75, las=3, col=4)
        axis(1, at=dom, labels=rep('',13))
        axis(2, at=((0:3)*0.2)/max_y1, labels=(0:3*0.2))

        # R0 panel
        plot(1:365, r0, col=orange, ylim=c(0,7), lwd=0.5, xlab='', axes=F, ylab=expression('R'[0]), type='n')
        shade_max=8
        shade_min=-1
#        polygon(c(season_start, season_end, season_end, season_start), c(shade_min, shade_min, shade_max, shade_max), col='#22222211',border=NA, xpd=T)
        polygon(c(1:365,365:1),c(r0_q975,rev(r0_q025)),col='#ffd6b2',border=NA)
#        polygon(c(season_start, season_end, season_end, season_start), c(shade_min, shade_min, shade_max, shade_max), col='#22222211',border=NA,xpd=T)
        abline(h=1, lty=2)
        legend('topright', legend=c(expression('Mean R'[0]),expression('Smoothed mean R'[0]),'95% IR'),
               col=c(1,orange,'#ffd6b2'), bty='n', lwd=c(0.5,2,0), pt.cex=2, lty=1, pch=c(NA,NA,15), cex=0.8)
#        lines(r0, lwd=1.5, col=orange)
        lines(r0_smooth, lwd=2, col=orange)
        lines(r0, lwd=0.5, col=1)
        box()
        axis(1, at=dom, labels=rep('',13))
        axis(2)
    }

    # EIP / mos / cases panel
    xlab_=''
    legend_loc = 'topleft'
    shade_max = 1.5
    ylim_ = c(0,1.5)
    if (supplement) {
        legend_loc = 'topright'
        shade_max = 1.11
    } else {
        ylim_ = c(0,1.0)
        par(mar=c(2.0, 2.1, 0, 1.1))
        par(oma=c(1.5,3,0.5,0.5))

        #par(mar=c(4.1, 4.1, 0.5, 1.1))
    }
    plot(rescaled_e_avg, ylim=ylim_, axes=F, ylab='', xlab=xlab_, type='n')
    #polygon(c(season_start, season_end, season_end, season_start), c(-0.135,-0.135,shade_max, shade_max), col='#22222222',border=NA,xpd=T)
    if(supplement) lines(r0_smooth/max(r0_smooth), col=orange, lwd=1.5)
    # EIP/temperature effect curve
    if(supplement) lines(rescaled_e_avg, col=3, lwd=1.5, type='l')

    if (supplement) {
        mtext("Month", 1, line=1, cex=0.75, outer=TRUE)
        #polygon(c(185,375,375,185), c(1.1,1.1,1.6,1.6), col='white', border=F)
        legend(legend_loc, legend=c('EIP', 'Mosquito population'), col=c(3,4), lty=1, lwd=1.5, bty='n', cex=0.8, inset=c(0.25,0))
        legend(legend_loc, legend=c(expression('Smoothed mean R'[0]),'Simulated cases','Observed cases'), col=c(orange,'#ee3333',1), lty=1, lwd=1.5, bty='n', cex=0.8)
        box()
    } else {
        box(bty='o')
        #text(195, 1.35, labels=expression(R[0]>=1))
        #legend(legend_loc, legend=c('EIP', 'Mosquito pop.','Simulated cases','Observed cases'), col=c(3,4,'#ee3333',1), lty=1, lwd=1.5, bty='n', cex=0.8, inset=c(0.015,0))
        legend(legend_loc, legend=c('Cases (observed)', 'Cases (model)'), col=c(1,'#ee3333'), lty=1, lwd=1.5, bty='n', cex=0.8, inset=c(0.07,0))
        text(0, 0.92, cex=2, font=2,labels='a')
    }
    # mosquitoes/precipitation effect curve
    if(supplement) {
        lines(1:365, seasonal_rain$mosquitoes/max(seasonal_rain$mosquitoes),col=4, lwd=1.5)
        axis(1, at=dom, labels=months)
    } else {
        axis(1, at=dom[seq(1,12,2)], labels=months[seq(1,12,2)], lwd = 0, line = -0.3)
        axis(1, at=dom, labels=F)
    }
    lines(week_days, epi, type='l', col=1, lty=1, lwd=1.5)
    lines(sim_epi/max(sim_epi), type='l', col='#ee3333', lwd=1.5)
    axis(2, at=0:2/2, labels=F, col=0, col.ticks=1)
    axis(2, at=0:3/2, labels=c(0.0, '', 'Max', ''), lwd.tick=0)
    #mtext("Scale", 2, line=0, cex=1, outer=TRUE)

    if (!supplement){
        mtext("Relative scale", side=2, outer=T, line=1, las=3)
        ylim_ = c(0,1.35)
        plot(rescaled_e_avg, ylim=ylim_, axes=F, ylab='', xlab='', type='n')
        box(bty='o')
        #lightgrey = '#ffffff'
        #darkgrey  = '#888888'
        lightgrey = '#eeeeff'
        darkgrey  = '#2222ff'
        gradient.rect(152,1.09,152+180,1.32,col=smoothColors(lightgrey,100,darkgrey,100,lightgrey),border=T)
        #gradient.rect(152,0,152+180,0.15,col=smoothColors(lightgrey,100,darkgrey,100,lightgrey),border=T)
        text(152+90, 1.20, labels='Insecticide active (optimal timing)')
        lines(rescaled_e_avg, col=3, lwd=1.5, type='l')
        lines(r0_smooth/max(r0_smooth), col=orange, lwd=1.5)
        abline(h=1/max(r0_smooth), col=orange, lty=3)
        lines(1:365, seasonal_rain$mosquitoes/max(seasonal_rain$mosquitoes),col=4, lwd=1.5)
        legend(legend_loc, legend=c('EIP', 'Mosquito pop.',expression('R'[0])), col=c(3,4,orange), lty=1, lwd=1.5, bty='n', cex=0.8, inset=c(0.07,0))
        text(0, ylim_[2]-(0.08*(ylim_[2]-ylim_[1])), cex=2, font=2,labels='b')

        axis(1, at=dom, labels=F)
        axis(1, at=dom[seq(1,12,2)], labels=months[seq(1,12,2)], lwd = 0, line = -0.3)
        axis(2, at=0:2/2, labels=F, col=0, col.ticks=1)
        axis(2, at=0:3/2, labels=c(0.0, '', 'Max', ''), lwd.tick=0)
        mtext("Month", 1, line=0, cex=1, outer=TRUE)
    }

    dev.off()
}

#plot_seasonality(supplement=T)
plot_seasonality(supplement=F)
