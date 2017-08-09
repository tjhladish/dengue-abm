require(data.table)
args <- commandArgs(trailingOnly = T)
# args <- c("~/Dropbox/who/fig1_data/foi_data.rds", "~/Dropbox/who/fig1_data/foi_effect.png")

eff <- readRDS(args[1])

foi = c('70% FOI', 'Baseline (fitted) FOI', '130% FOI')
seasons = c('June 1 IRS (best timing)', 'November 1 IRS (worst timing)')

data = matrix(setkey(eff, timing)$med.eff, nrow=2, byrow=T)
rownames(data) = seasons
colnames(data) = foi
png(args[2], width=1200, height=800, res=120)
barplot(data, beside=T, ylab='Overall cumulative effectiveness (first 10 years)', legend=rownames(data), main='Effect of force of infection and IRS seasonality on effectiveness', ylim=c(0,1))
dev.off()
