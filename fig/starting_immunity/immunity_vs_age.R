for (i in 1:12) {
    d = read.table(paste0('imm_vs_age.', i))

    png(paste0('age_vs_immunity', i, '.png'), height=800, width=1000, res=150)
    par(mar=c(5,4.5,2,2))
    plot(d$V1, d$V3/d$V2, type='l', ylim=c(0,1), xlab='Age (years)', ylab='Fraction immune', lwd=2, col='blue')
    points(d$V1, d$V4/d$V2, type='l', ylim=c(0,1), col='red', lwd=2)
    points(d$V1, d$V5/d$V2, type='l', ylim=c(0,1), col='green', lwd=2)
    points(d$V1, d$V6/d$V2, type='l', ylim=c(0,1), col='purple', lwd=2)
    legend('top',fill=c('blue','red','green','purple'),legend=c('DENV1     ','DENV2     ','DENV3     ','DENV4'), horiz=T, bty='n')
    dev.off()
}
