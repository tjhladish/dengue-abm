require(data.table)
args <- commandArgs(trailingOnly = T)
# args <- c("~/Dropbox/who/fig1_data/foi-eff.rds","~/Dropbox/who/fig1_data/foi_effect.png")

eff <- readRDS(args[1])

#foi = 
seasons = c('June 1 IRS (proactive timing)', 'November 1 IRS (reactive timing)')

data = matrix(setkey(eff, doy)[year==9,cq.med], nrow=2, byrow=T)
rownames(data) = seasons

png(tail(args,1), width=1200, height=800, res=120)
barplot(data,
  beside=T,
  ylab='Overall cumulative effectiveness (first 10 years)',
  names.arg = c(expression(paste('70% ',M[peak])), expression(paste('Baseline (fitted) ',M[peak])), expression(paste('130% ',M[peak]))),
  legend=rownames(data),
  main='Effect of mosquito population and IRS seasonality on effectiveness', ylim=c(0,1)
)
dev.off()
