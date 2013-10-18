d = read.table('yucatan.edgelist', header=TRUE)
f = d[1:100000,]
plot(c(f$x1, f$x2), c(f$y1, f$y2), pch='.')
segments(f$x1, f$y1, f$x2, f$y2, col='blue')

library('sp')
m = as.matrix(f)

# pixel width in km
spDists(matrix(c(-88.29992664, 20.63328984),byrow=1,nrow=1), matrix(c(-88.29575997, 20.63328984),byrow=1,nrow=1), longlat=TRUE)

# pixel height in km
spDists(matrix(c(-88.28325996, 20.63328984),byrow=1,nrow=1), matrix(c(-88.28325996, 20.62912317),byrow=1,nrow=1), longlat=TRUE)

m=as.matrix(d)
distances = apply(m, 1, function(x) spDists(matrix(x[1:2],byrow=1,nrow=1), matrix(x[3:4],byrow=1,nrow=1), longlat=TRUE))
median(distances)
mean(distances)
min(distances)
hist(distances[distances<.5],nclass=100)
