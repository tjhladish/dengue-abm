require(data.table)
obs = read.table('weekly_confirmed_cases_1995-2015_wide.csv', header=F, sep=' ', col.names=1:53, fill=T)
#obs = obs[-c(1:11),]                   # discard 1995-2015
# sim.raw = read.table('daily_output-refit', header=F, col.names=c('serial','day','intro','local','intro_prev','local_prev'))
# sim.raw = sim.raw[sim.raw$day<=50004,] # end before 2016
# #sim.raw = sim.raw[sim.raw$day>46354,]  # start with 2006
# sim.raw$inc = sim.raw$intro + sim.raw$local
# sim.dt = data.table(sim.raw)[,.(serial,day,inc)]
# sim.means = aggregate(sim.dt, by = list(serial=sim.dt$serial,day=sim.dt$day%%365), FUN=mean)
#
# sim.dt.mean = data.table(sim.means)[,.(serial,day,inc)]
load(file='sim.dt.mean.rda')

sim.overall.mean = aggregate(sim.dt.mean[,.(inc)], by=list(day=sim.dt.mean$day), FUN=mean)
obs.mean = colMeans(obs[,1:52], na.rm=T)

sim.ci = aggregate(sim.dt.mean[,.(inc)], by=list(day=sim.dt.mean$day), FUN=function(x) quantile(x, probs=c(0.05, 0.95)))

png('dengue_seasonality_1995-2015-refit.png', height=600, width=1200, res=120)
#plot(seq(1,365,7), obs.mean[c(1:52,1)], type='l', ylim=c(0,1.5), main='Dengue seasonality in Yucatan (2006-2015)', ylab='Incidence', xlab='Julian day',lwd=2)
plot(seq(1,365,7), obs.mean[c(1:52,1)]/max(obs.mean), type='l', ylim=c(0,1.25), main='Dengue seasonality in Yucatan (1995-2015)', ylab='Incidence', xlab='Julian day',lwd=3)
#matlines((1:53)*7, t(obs), type='l', lty=1, lwd=0.5)
matlines(seq(1,363,7), t(obs[15:19,1:52]/obs[cbind(15:19,max.col(obs[15:19,1:52]))]), type='l', lwd=0.5, lty=1, col=c('red','orange','green','purple','darkgrey'))
legend('topleft',legend=c('Observed','Simulated'),col=c('black','blue'), lwd=2, bty='n')
lines((1:365), sim.overall.mean$inc/max(sim.overall.mean$inc), col='blue', lwd=3)
polygon(c(1:365, 365:1), c(sim.ci$inc[,1],rev(sim.ci$inc[,2]))/max(sim.overall.mean$inc), col = "#00009922", border = FALSE)
dev.off()

