eff = c(0.8160120605, 0.6574398432, 0.5511642524, 0.5837718242, 0.497082237, 0.3995614269)
foi = c('70% FOI', 'Baseline (fitted) FOI', '130% FOI')
seasons = c('June 1 IRS (best timing)', 'November 1 IRS (worst timing)')
data = table(eff[1:3], eff[4:6])

data = matrix(eff, nrow=2, byrow=T)
rownames(data) = seasons
colnames(data) = foi
png('foi_effect.png', width=1200, height=800, res=120)
barplot(data, beside=T, ylab='Overall cumulative effectiveness (first 10 years)', legend=rownames(data), main='Effect of force of infection and IRS seasonality on effectiveness', ylim=c(0,1))
dev.off()
