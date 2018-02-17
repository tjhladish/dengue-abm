rm(list=ls())

require(data.table)
require(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
# args <- c(paste0("~/Dropbox/who/fig1_data/stopping-","sero-SI",".rds"),"~/Dropbox/who/fig1_data/SI-fig-sero.png")

# three cases
#  no intervention
#  75% coverage
#  75% coverage w/ stopping
seroprevalence.dt <- readRDS(args[1])

idvars = c("coverage", "duration", "end_year", "year")

sero.means <- seroprevalence.dt[,
  .SD,
  by=idvars,
  .SDcols=grep("\\.mn", names(seroprevalence.dt), value = T)
]
# [,
#   imm1.mn := imm1.mn + imm0.mn,
#   by=idvars
# ][,
#   imm2.mn := imm2.mn + imm1.mn,
#   by=idvars
# ][,
#   imm3.mn := imm3.mn + imm2.mn,
#   by=idvars
# ][,
#   imm4.mn := 1,
#   by=idvars
# ]

sero.mlt <- melt(sero.means, id.vars = idvars)

sero.mlt[, prev_infs := as.integer(gsub("imm(\\d).+","\\1",variable))]

sero.baseline <- setkey(rbind(
  data.table::copy(sero.mlt[end_year == 0])[, coverage := 25 ],
  data.table::copy(sero.mlt[end_year == 0])[, coverage := 50 ],
  data.table::copy(sero.mlt[end_year == 0])[, coverage := 75 ]
), prev_infs
)

ggplot(
  rbind(sero.baseline, sero.mlt[end_year!=0]),
  aes(
    x=year, y=value,
    linetype=factor(end_year), color=factor(prev_infs),
    group=interaction(variable, end_year)
  )
) + geom_line() + theme_minimal() + facet_grid(coverage ~ .)

# plot_years = 20  # used 5 for timing plot
# 
# .col = c('#cce3ff','#99c7ff','#4a91cc','#2f5b7f','#104060')
# stackimmunity = function(r,.label='') {
#     .range = 11:30
#     mat = matrix(c(mean0[r,.range],mean1[r,.range],mean2[r,.range],mean3[r,.range],mean4[r,.range]), ncol=5, nrow=length(.range))
#     stackpoly(mat, stack=T, ylim=c(0,1), col=.col, axis4=F, xaxlab = F, axis=F)
#     axis(2,tick = T,labels=c('0.0','0.2','0.4','0.6','0.8','1.0'),at=0:5/5, las=1)
#     axis(4,labels = F)
#     text(1.2, 0.07, pos=4, labels=.label, font=2)
# }
# 
# #legend('topright',legend = c('all houses treated in 1 day', 'houses treated across 90 days', 'houses treated continuously'), lty=c(2,1,3), lwd=lwd_, bty='n')
# 
# png('susceptible_and_mean_immunity.png',height=1200,width=1000,res=180)
# #mean0 = aggregate(data0, by=list(tags$campaign_duration, tags$vc_coverage, tags$timing, tags$strat_years, tags$vector_control), FUN=mean)
# rows = c(1,7,4)
# .range = 11:30
# par(mfrow=c(2,1),mar=c(1,4.2,2,1),oma=c(3.5,1,0,0))
# matplot(t(mean0[rows,.range]), type='l', col='black', ylim=c(0,0.6), ylab='Fully susceptible population', xlab='')
# legend('bottomright',inset=c(0.025,0.05),legend=c('Baseline (No IRS)','IRS, 75% coverage','IRS, 75% coverage ended in year 10'), lty=1:3,bty='n',cex=0.75)
# avg_inf = t(mean1[rows,.range] + 2*mean2[rows,.range] + 3*mean3[rows,.range] + 4*mean4[rows,.range])
# matplot(avg_inf, type='l', col='black', ylim=c(0,1.5), ylab='Mean number of past infections', xlab='')
# mtext("Intervention year", side=1, outer=T, line=1.5)
# dev.off()
