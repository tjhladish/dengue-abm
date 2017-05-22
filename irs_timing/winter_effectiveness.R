rm(list=ls())
require("RSQLite")
require("data.table")
drv = dbDriver("SQLite")
path <- sprintf("%s%s",.dbox.path,"who")
db.path <- sprintf("%s/%s",path,"irs_timing-summer_winter-foi.sqlite")
log.path <- sprintf("%s/%s",path,"daily-irs_refit-simple-nofilenames-uniq.out")

db = dbConnect(drv, db.path)

# campaign_duration is always 1
# coverage is always 75%

ref.dt <- data.table(
  dbGetQuery(db,
    'SELECT P.serial AS serial,
       vector_control,
       CASE timing WHEN 152 THEN "summer" ELSE "winter" END AS intervention,
       CASE WHEN foi=1.0 THEN "fit" WHEN foi<1.0 THEN "low" ELSE "high" END AS foi
     FROM par P, job J
       WHERE P.serial = J.serial
         AND status = \'D\';'
  )
)

ref.dt[, foi := factor(foi, levels = c("low","fit","high"), ordered = T)]
ref.dt[vector_control == 0, intervention := "none" ]
ref.dt$vector_control <- NULL

daily.dt <- fread(log.path,sep = " ", header = F)
setnames(daily.dt, c("V1","V2","V3","V4","V5","V6","V7"),
         c("V1","serial","day","incidence.intro","prevalence.intro","incidence.local","prevalence.local"))
daily.dt$V1 <- NULL

require(reshape2)

plot.dt <- setkey(melt.data.table(daily.dt, id.vars=c("serial","day")),serial,day)[ref.dt, on="serial"][,{
  mn <- mean(value)
  ps <- quantile(value, probs = c(0,.25,.5,.75, 1))
  list(ave=mn, min=ps[1], lo=ps[2], md=ps[3], hi=ps[4], max=ps[5])
}, keyby=list(intervention, foi, variable, day)
]
rm(daily.dt)
plot.dt[, day    := day - min(day) ]


require(ggplot2)

ggplot(plot.dt[variable == "prevalence.local"],
  aes(x=day,color=intervention,fill=intervention)
) + theme_minimal() +
  coord_cartesian(xlim=c(125,575), ylim=c(1,10000)) +
  facet_grid(. ~ foi, scales = "free") +
#  geom_ribbon(aes(ymin=min, ymax=max, alpha="full range", color=NULL)) +
  geom_ribbon(aes(ymin=lo, ymax=hi, alpha="half range", color=NULL)) +
#  geom_line(aes(y=ave, linetype="average", alpha="line")) +
  geom_line(aes(y=md, linetype="median", alpha="line")) +
  scale_y_log10(name="absolute prevalence (log scaling)")

npars = 4
serial_col = npars + 1
data_burnin = 5 # used 6 for timing plot
plot_years = 25  # used 5 for timing plot
last_col = serial_col + data_burnin + plot_years

tags = d[,1:npars]
data = d[,(serial_col + data_burnin + 1):last_col]
#medians = aggregate(d[,5:34], by=list(d$vector_control, d$vc_coverage, d$timing), FUN=median)

plot_effectiveness_over_time = function(tags, data, timing) {
  
  means = aggregate(data, by=list(tags$vector_control, tags$campaign_duration, tags$vc_coverage, tags$timing), FUN=median)
  names(means)[1:npars]=c('vc','dur','cov','day')
  means = means[means$day==timing | means$vc==0,]
  eff_vals = 1 - sweep(as.matrix(means[-1,(npars+1):(dim(means)[2])]), 2, as.numeric(means[1,(npars+1):(dim(means)[2])]), '/')
  eff = cbind(means[-1,1:npars], eff_vals)
  
  cum_means = cbind(means[,1:npars], t(apply(means[,-c(1:npars)], 1, cumsum)))
  cum_eff_vals = 1 - sweep(as.matrix(cum_means[-1,(npars+1):(dim(cum_means)[2])]), 2, as.numeric(cum_means[1,(npars+1):(dim(cum_means)[2])]), '/')
  cum_eff = cbind(cum_means[-1,1:npars], cum_eff_vals)
  
  #par(mfrow=c(2,1), mar=c(2.1,2.1,1,1), oma=c(3,3,3,0))
  matplot(t(eff[eff$dur==1,5:29]), lty=3:1, type='l', ylim=c(-0.2,1), ylab='Effectiveness',col='#DD0000')
  abline(h=0,lty=2)
  #matplot(t(cum_eff[cum_eff$dur==1,5:29]), lty=3:1, type='l', ylim=c(0,1), ylab='Cumulative effectiveness',lwd=2,main='Cumulative overall effectiveness of IRS',xlab='Year',col='#DD0000')
  legend('topright',legend=c('75% coverage','50% coverage','25% coverage'),lwd=2,lty=1:3,bty='n',col='#DD0000')
  return(means)
}

png('vc_effectiveness.png', width=1200, height=960, res=150)
par(las=1,bty='L')
plot_effectiveness_over_time(tags,data,182)
dev.off()

timing_effectiveness = function(tags, data, window, func, end) {
  if (end) {
    d = cbind(tags, rowSums(data[,(dim(data)[2]- window + 1):dim(data)[2]]))
  } else {
    d = cbind(tags, rowSums(data[,1:window]))
  }
  means = aggregate(d[,dim(d)[2]], by=list(d$vector_control, d$campaign_duration, d$vc_coverage, d$timing), FUN=func)
  names(means)=c('vc','dur','cov','day','sum')
  means$eff = 1 - means$sum/means$sum[1]
  return(means)
}

png('ce_vc-equillibrium-90day-first5years.png', width=12, height=6.4, res=200, units='in')
par(mfrow=c(1,2), mar=c(2.1,2.1,1,1), oma=c(3,3,3,0),las=1)
#dur_ = 1
end = F
for (cov_ in c(0.25, 0.75)){
  test = timing_effectiveness(tags, data, plot_years, mean, end)
  plot(0:365, ylim=c(0,1), xlab='', ylab='', type='n')
  minval = min(c(test$eff[test$vc==1 & test$dur==0 & test$cov==cov_],test$eff[test$vc==1 & test$dur==1 & test$cov==cov_ & test$day < 365]))
  maxval = max(c(test$eff[test$vc==1 & test$dur==0 & test$cov==cov_],test$eff[test$vc==1 & test$dur==1 & test$cov==cov_ & test$day < 365]))
  rect(0,minval,365,maxval,border=NA,col='#cccccc')
  lines(test$day[test$vc==1 & test$dur==0 & test$cov==cov_], test$eff[test$vc==1 & test$dur==0 & test$cov==cov_], ylim=c(0,1), type='l', xlab='', ylab='', lty=2)
  
  #    x90 = ((45+test$day[test$vc==1 & test$dur==1 & test$cov==cov_ & test$day < 365])%%365)[c(47:53,1:46)]
  x90 = (test$day[test$vc==1 & test$dur==1 & test$cov==cov_ & test$day < 365])
  y90 = (test$eff[test$vc==1 & test$dur==1 & test$cov==cov_ & test$day < 365])
  lines(x90, y90)
  legend('bottomright',inset = c(0.08,0.1), legend = paste0(cov_*100, '% household coverage'), bty='n')
  legend('bottomright',legend = c('all houses treated in 1 day', 'houses treated across 90 days'), lty=c(2,1), bty='n')
}
par(las=0)
mtext("Campaign start (Julian day)", side=1, outer=T, line=1)
mtext("Effectiveness (prevented infections)", side=2, outer=T, line=1)
mtext("Effect of campaign timing on IRS effectiveness", side=3, outer=T)
dev.off()
