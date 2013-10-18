#d = read.table('yucatan_edges_w_lengths.out',header=T)
#a = read.table("../../raw_data/all_yucatan_pixels.out", sep='\t')


f = d[d$dist_km <= 10.0,]
#f = f[1:100000,]

makeplot = function() {
    png('edges.png',width=1200, height=1000, res=130)
    par(mar=c(5,5,2,2))
    plot(a$V1, a$V2, pch='.', ylim=c(19.6, 21.65), xlim=c(-90.4, -87.5), asp=1, xlab='Longitude', ylab='Latitude', col='black', cex=1.5, cex.lab=1.5)

    #plot(c(f$x1, f$x2), c(f$y1, f$y2), pch='.', asp=1, cex=1.5, cex.lab=1.5, xlab='Longitude', ylab='Latitude',)
    segments(f$x1, f$y1, f$x2, f$y2, col='red')
}
