monthly_ts = read.table('imm_seq_100_years_hfoi/monthly_ts_by_rank')

for (r in c(1)) {
#for (r in c(1,2,3,4,7,8,9,11)) {
#    pdf(paste0('age_vs_immunity_sequence.', r, '.pdf'), height=8, width=10)
    idx = 0
    for (y in 0:99) {
        for (m in 0:12) {
            if (m == 0 && y > 0) next;
            idx = idx + 1
            d = read.table(paste0('imm_seq_100_years_hfoi/immunity.', r, '.year', y, '.mon', m))
            png(paste0('age_vs_immunity.', r, '.', idx, '.png'), width=1200, height=900, res=120)

            #png(paste0('age_vs_immunity', i, '.png'), height=1000, width=1000, res=150)
            layout(matrix(c(1,2)), heights=c(1,2))
            par(mar=c(3,4.5,2,2))
            ts_data = as.numeric(monthly_ts[monthly_ts$V1==1,-1:-12])
            plot(ts_data, type='l', main='', xlab='', ylab='Incidence', axes=F)
            mtext('Year', side=1, line=2)
            axis(1, at=seq(1,1201,120), labels=seq(0,100,10))
            axis(2)
            points(idx, ts_data[idx], col='#FF000099', cex=2, pch=20)
            par(mar=c(5,4.5,1,2))
            plot(d$V1, d$V3/d$V2, type='l', ylim=c(0,1), col='red', lwd=2,
                 xlab='', ylab='Fraction immune', axes='F')
            mtext('Age', side=1, line=2)
            axis(1); axis(2)
            points(d$V1, d$V4/d$V2, type='l', ylim=c(0,1), col='green', lwd=2)
            points(d$V1, d$V5/d$V2, type='l', ylim=c(0,1), col='blue', lwd=2)
            points(d$V1, d$V6/d$V2, type='l', ylim=c(0,1), col='purple', lwd=2)
            points(d$V1, d$V7/d$V2, type='l', ylim=c(0,1), col='grey', lwd=2)
            points(d$V1, d$V8/d$V2, type='l', ylim=c(0,1), col='black', lwd=2)
            legend('topleft',
                   fill=c('red','green','blue','purple','grey','black'),
                   legend=c('DENV1','DENV2','DENV3','DENV4', 'Any', 'All'), bty='n', inset=0.02)
            dev.off()
        }
    }
#    dev.off()
}

