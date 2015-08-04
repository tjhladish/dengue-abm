rm(list=ls())
d = read.table('delaunay_edge_data.out',header=T)
a = read.table("../../raw_data/all_yucatan_pixels.out", sep='\t')
library(plotrix)

f = d[d$d_km <= 1,]
tbox_ht = 0.0045

makeplots = function() {
    png('edges_1km.png',width=1200, height=1000, res=130)
    par(mar=c(5,5,2,2))
    plot(a$V1, a$V2, pch='.', ylim=c(19.6, 21.65), xlim=c(-90.4, -87.5), asp=1, xlab='Longitude', ylab='Latitude', col='#dddddd', cex=1.5, cex.lab=1.5)
    rect(-89.035, 20.94, -88.905, 21.04, col="#ffffff", border = NA)

    segments(f$x1, f$y1, f$x2, f$y2, col='red')
    rect(-89.035, 20.94, -88.905, 21.04, col="#ffffff00")
    lines(c(-88.0, -87.5229), c(19.6, 19.6), lwd=3)
    text(-87.75, 19.55, "50 km")
    dev.off()

    png('edges_1km-detail.png',width=1200, height=1000, res=130)
    par(mar=c(5,5,2,2))
    plot(a$V1, a$V2, pch='.', ylim=c(20.945,21.035), xlim=c(-89.02, -88.92), asp=1, xlab='Longitude', ylab='Latitude', col='black', cex=1.5, cex.lab=1.5, type='n')
    segments(f$x1, f$y1, f$x2, f$y2, col='red')

    textbox(c(-88.935,-88.935), 21.0325, "Tekal de Venegas", margin=c(0,0,0.001,0.0005))

    textbox(c(-88.946,-88.946), 20.958, "Sitilpech", margin=c(0,0,0.001,0.0005))

    rect(-89.027, 20.9455, -89.016, 20.9455+tbox_ht, col="#ffffff")
    textbox(c(-89.027, -89.016), 20.95, "Izamal", margin=c(0,0,0.001,0),box=F)

    lines(c(-88.92963, -88.92), c(20.9505, 20.9505), lwd=3)
    text(-88.925, 20.9485, "1 km")
    dev.off()
}

makeplots()
