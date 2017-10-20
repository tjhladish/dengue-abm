rm(list=ls())

require(data.table)
require("plotrix")

args <- commandArgs(trailingOnly = T)
# args <- c("~/Dropbox/who/fig1_data/stopping-baseline.rds", "~/Dropbox/who/fig1_data/stopping-interventions.rds", "~/Dropbox/who/fig1_data/stacked_immunity.png")

baseline.sero <- readRDS(args[1])
interventions.sero <- readRDS(args[2])

stat.base <- baseline.sero[,.(
  mn.0=mean(imm0),
  mn.1=mean(imm1),
  mn.2=mean(imm2),
  mn.3=mean(imm3),
  mn.4=mean(imm4),
  coverage = 0,
  end_year = 50
), by=year
]

stat.int <- interventions.sero[,.(
  mn.0=mean(imm0),
  mn.1=mean(imm1),
  mn.2=mean(imm2),
  mn.3=mean(imm3),
  mn.4=mean(imm4)
), by=.(coverage, end_year, year)
]

dt <- rbind(stat.int, stat.base)

.col = c('#cce3ff','#99c7ff','#4a91cc','#2f5b7f','#104060')
stackimmunity = function(dt, scen, maxyrs, .label='') {
  slice <- dt[eval(scen) & year < maxyrs,.SD,.SDcols=paste0("mn.",0:4)]
    mat = as.matrix(slice)
    stackpoly(mat, stack=T, ylim=c(0,1), col=.col, axis4=F, xaxlab = F, axis=F)
    axis(2,tick = T,labels=c('0.0','0.2','0.4','0.6','0.8','1.0'),at=0:5/5, las=1)
    axis(4,labels = F)
    text(1.2, 0.07, pos=4, labels=.label, font=2)
}

#legend('topright',legend = c('all houses treated in 1 day', 'houses treated across 90 days', 'houses treated continuously'), lty=c(2,1,3), lwd=lwd_, bty='n')

png(args[3],height=1200,width=1000,res=180)
par(mfrow=c(3,1),mar=c(1.5,3,0.5,1),oma=c(3,2.5,0.5,0.5),xpd=T)
stackimmunity(dt,expression(coverage==0),20,'Baseline (No IRS)')
legend('bottomright',inset=c(0.025,0.05),legend=c('4  ','3  ','2  ','1  ','0  '),title='Number of past infections',fill=rev(.col),horiz = T,x.intersp = 0.75)
stackimmunity(dt,expression(coverage==75 & end_year==50),20,'IRS, 75% coverage')
stackimmunity(dt,expression(coverage==75 & end_year==10),20,'IRS, 75% coverage ended in year 10')
axis(1,tick = F,labels=1:4*5,at=1:4*5)

mtext("Intervention year", side=1, outer=T, line=1.5)
mtext("Population immunity", side=2, outer=T, line=0.5)
#mtext("Effect of intervention on population immunity", side=3, outer=T)
dev.off()