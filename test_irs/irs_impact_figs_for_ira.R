rm(list=ls())

require("RSQLite")
drv = dbDriver("SQLite")
db = dbConnect(drv, "./vacvec-10x.sqlite")
basic <- dbGetQuery(db, 'select vector_control, vc_coverage, catchup, vac, catchup_to, M.*
                      from parameters P, metrics M, jobs J
                      where P.serial = M.serial 
                            and P.serial = J.serial
                            and foi_target = 2
                            and target = 9
                            and vc_coverage = 0.25
                            and status = \'D\';')

tag_cols = 1:5
tags = basic[,tag_cols]
data = basic[,12:36]

medians = aggregate(data, by=list(tags$vac, tags$catchup, tags$catchup_to, tags$vector_control, tags$vc_coverage), FUN=median)
names(medians)[tag_cols] = c("vac", "catchup", "catchup_to", "vector_control", "vc_coverage")

#labels = d[,1:3]
#vals = d[,4:28]
cum_data = t(apply(data, 1, cumsum))
names(cum_data) = names(data)
cd = cbind(tags, cum_data)
cum_medians = aggregate(cd, by=list(cd$vac, cd$catchup, cd$catchup_to, cd$vector_control, cd$vc_coverage), FUN=median)[-tag_cols]

png('dengue_cases-rerun10x.png', width=9, height=7, res=200, units='in')
matplot(t(medians[medians$vector_control==0,-tag_cols]/8.39660), type='l', lty=c(1,1,3,1), lwd=2, ylab='Annual incidence (cases per 100k)', col=c('black', 'green3',4,4), xlab='Year', main='Impact of vaccination on dengue cases (WHO assumptions -- 80% coverage)')
legend('bottomright',bty='n',legend = c('no vaccine','routine 9-year-old','routine 9yo + catchup to 17','routine 9yo + catchup to 30'), col=c(1,'green3',4,4), lty=c(1,1,3,1), lwd=2)
dev.off()

effect = 1 - sweep(as.matrix(medians[medians$vac==1 & medians$vector_control==0,-tag_cols]), 2, as.numeric(medians[1,-tag_cols]), '/' )

png('dengue_effectiveness-cases-rerun10x.png', width=9, height=7, res=200, units='in')
matplot(t(effect), type='l', lty=c(1,3,1), lwd=2, ylab='Annual effectiveness', col=c('green3',4,4), xlab='Year', main='Effectiveness of vaccination (WHO assumptions -- 80% coverage)')
legend('bottomright',bty='n',legend = c('routine 9-year-old','routine 9yo + catchup to 17','routine 9yo + catchup to 30'), col=c('green3',4,4), lty=c(1,3,1), lwd=2)
abline(h=0, lty=3, lwd=0.75)
dev.off()


png('dengue_cumul_cases-rerun10x.png', width=9, height=7, res=200, units='in')
matplot(t(cum_medians[cum_medians$vector_control==0,-tag_cols]/8.39660), type='l', lty=c(1,1,3,1), lwd=2, ylab='Cumulative incidence (cases per 100k)', col=c('black', 'green3',4,4), xlab='Year', main='Cumulative impact of vaccination (WHO assumptions -- 80% coverage)')
legend('bottomright',bty='n',legend = c('no vaccine','routine 9-year-old','routine 9yo + catchup to 17','routine 9yo + catchup to 30'), col=c(1,'green3',4,4), lty=c(1,1,3,1), lwd=2)
dev.off()

cum_effect = 1 - sweep(as.matrix(cum_medians[-1,-tag_cols]), 2, as.numeric(cum_medians[1,-tag_cols]), '/' )

png('dengue_cumul_effectiveness-VC-routine_only-cases-rerun10x.png', width=9, height=7, res=200, units='in')
matplot(t(cum_effect[c(1,4,5),]), type='l', lty=1, lwd=c(2,2,3), ylab='Cumulative effectiveness', col=c(4,2,'purple'), xlab='Year', main='Cumulative effectiveness of vector control and routine vaccination ', ylim=c(-0.25,1))
legend('topright',bty='n',legend = c('routine 9-year-old','routine IRS (vector control)','both'), col=c(4,2,'purple'), lty=1, lwd=c(2,2,3))
abline(h=0, lty=3, lwd=0.75)
dev.off()

png('dengue_cumul_effectiveness-VC-w_catchup-cases-rerun10x.png', width=9, height=7, res=200, units='in')
matplot(t(cum_effect[c(3,4,7),]), type='l', lty=1, lwd=c(2,2,3), ylab='Cumulative effectiveness', col=c(4,2,'purple'), xlab='Year', main='Cumulative effectiveness of vector control and routine vaccination + catchup', ylim=c(-0.25,1))
legend('topright',bty='n',legend = c('routine 9yo + catchup to 30', 'routine IRS (vector control)','both'), col=c(4,2,'purple'), lty=1, lwd=c(2,2,3))
abline(h=0, lty=3, lwd=0.75)
dev.off()

ca = sweep(-as.matrix(medians[medians$vac==1 & medians$vector_control==0,-tag_cols]), 2, as.numeric(medians[1,-tag_cols]), '+' )
VAR = ca/839660
