rm(list=ls())
require(data.table)
require(raster)
require(maptools)
require(magrittr)

mex = shapefile("/home/tjhladish/DropboxUFL-data/areas_geoestadisticas_estatales.shp") %>%
      spTransform(CRS("+init=epsg:4326"))
yuc = shapefile("/home/tjhladish/DropboxUFL-data/shapefiles_INEGI_SCINCE_2010/estatal.shp") %>%
      spTransform(CRS("+init=epsg:4326"))
mun = shapefile("/home/tjhladish/DropboxUFL-data/shapefiles_INEGI_SCINCE_2010/municipal.shp") %>%
      spTransform(CRS("+init=epsg:4326"))

# Get just the mainland polygon for Yucatan
yuc <- yuc@polygons[[1]]@Polygons[[7]] %>%
    list() %>%
    Polygons(ID=1) %>%
    list() %>%
    SpatialPolygons()
    #SpatialPolygons(proj4string = CRS("+init=epsg:4326"))
#plot(ntl_ras)

pop = fread('/home/tjhladish/work/dengue/pop-yucatan/population-yucatan.txt')
loc = fread('/home/tjhladish/work/dengue/pop-yucatan/locations-yucatan.txt')
net = fread('/home/tjhladish/work/dengue/pop-yucatan/network-yucatan.txt')

all_houses  = subset(loc, type == 'h')
all_schools = subset(loc, type == 's')
all_works   = subset(loc, type == 'w')

mos_movements = merge(net, loc, by.x='locationID1', by.y='id')
mos_movements = merge(mos_movements, loc, by.x='locationID1.1', by.y='id')
names(mos_movements) = c('locid1', 'locid2','type1','x1', 'y1', 'type2', 'x2', 'y2')


#plot(yuc, axes=T, col='black') # Can we plot a bassic shape?
#attach(mos_movements)
#segments(x1, y1, x2, y2, col='green')
#detach(mos_movements)
#points(all_houses$x, all_houses$y, pch='.', col='#b23f2b')
#points(all_houses$x, all_houses$y, pch='.', col='black')
#points(all_works$x, all_works$y, pch=20, col='#216ea6', cex=0.75)
#points(all_schools$x, all_schools$y, pch=20, col='#dd0000')


.xall = c(-90.52197, -87.41819)
.yall = c(19.46822,  21.70804)
#.xmed = c(-89.52,-89.45)
#.ymed = c(20.70,20.78)
#.xlrg = c(-89.44,-89.31) # Tonkal?
#.ylrg = c(20.565,20.895)
#.xmed = c(-89.49,-89.46)
#.ymed = c(20.715,20.745)

.med_width = 0.01
 
#ctr = c(20.965638, -88.603060) # cenotillo
ctr = c(20.932985, -89.018394)  # izamal

.xmed = c(ctr[2] - .med_width, ctr[2] + .med_width)
.ymed = c(ctr[1] - .med_width, ctr[1] + .med_width)

#.xmed = c(-89.626, -89.616) #old
#.ymed = c(20.882, 20.89)    #old
.sml_width = 0.001

#ctr = c(20.933392, -89.025688)
ctr = c(20.930862, -89.024314)

.xsml = c(ctr[2] - .sml_width, ctr[2] + .sml_width)
.ysml = c(ctr[1] - .sml_width, ctr[1] + .sml_width)

.xlim = .xall
.ylim = .yall

.xlim = .xmed
.ylim = .ymed

.xlim = .xsml
.ylim = .ysml



draw_mosquito_movements = function(mos_movements, .xlim, .ylim) {
    .epsilon = 0.02 # aprox 2 km in degrees
    mos_movements_subset = subset(mos_movements, x1 >= .xlim[1] - .epsilon &
                                      x2 <= .xlim[2] + .epsilon &
                                      y1 >= .ylim[1] - .epsilon &
                                      y2 <= .ylim[2] + .epsilon)

    attach(mos_movements_subset)
    #segments(x1, y1, x2, y2, col='#90dd80')
    segments(x1, y1, x2, y2, col='#619e6a80')
    detach(mos_movements_subset)
}

get_people_movements = function(.xlim, .ylim) {
    houses  = subset(loc, x >= .xlim[1] & x <= .xlim[2] & y >= .ylim[1] & y <= .ylim[2] & type == 'h')
    schools = subset(loc, x >= .xlim[1] & x <= .xlim[2] & y >= .ylim[1] & y <= .ylim[2] & type == 's')
    works   = subset(loc, x >= .xlim[1] & x <= .xlim[2] & y >= .ylim[1] & y <= .ylim[2] & type == 'w')
    
    pop_subset = subset(pop, hid %in% houses$id | workid %in% c(schools$id, works$id))
    pop_movements = subset(pop_subset, hid != workid)
    pop_movements = merge(pop_movements, loc, by.x='hid', by.y='id')
    pop_movements = merge(pop_movements, loc, by.x='workid', by.y='id')
    pop_movements = pop_movements[,c('pid', 'hid', 'age', 'sex', 'workid', 'empstat', 'type.y', 'x.x', 'y.x', 'x.y', 'y.y')]
    pm_col_ct = length(pop_movements)
    names(pop_movements)[pm_col_ct - 4:0] = c('type.work', 'x1', 'y1', 'x2', 'y2')
    
    # plot just houses, workplaces, & schools within .xlim and .ylim
    # points(houses$x, houses$y, pch='.', col='black')
    # points(works$x, works$y, pch=20, col='blue', cex=0.5)
    # points(schools$x, schools$y, pch=20, col='red')

    return(pop_movements)
}

make_village_map = function(make_png=F, ctr = c(20.932985, -89.018394), .med_width = 0.01) {
    .xlim = c(ctr[2] - .med_width, ctr[2] + .med_width)
    .ylim = c(ctr[1] - .med_width, ctr[1] + .med_width)
    
    if (make_png) {
        png('izamal_model-recolor-legend.png', width=2400, height=1600, res=240)
        par(mar=c(5.1,4.1,2,2))
    }
    
    plot(yuc, axes=T, xlim=.xlim, ylim=.ylim, xlab='Longitude', ylab='Latitude', col='#e3debf80')
    draw_mosquito_movements(mos_movements, .xlim, .ylim)

    pop_movements = get_people_movements(.xlim, .ylim)    

    student_movements = subset(pop_movements, type.work == 's')
    work_movements = subset(pop_movements, type.work == 'w')
    
    # plot how workers and students move?
    #attach(work_movements)
    ##rand_choice = runif(length(x1)) < 0.01
    #rand_choice = T
    #segments(x1[rand_choice], y1[rand_choice], x2[rand_choice], y2[rand_choice], col='#000000')
    #detach(work_movements)
    
    #attach(student_movements)
    ##rand_choice = runif(length(x1)) < 0.01
    #rand_choice = T
    #segments(x1[rand_choice], y1[rand_choice], x2[rand_choice], y2[rand_choice], col='#ff0000')
    #detach(student_movements)

    points(all_houses$x, all_houses$y, pch='.', col='#b23f2b')
    points(all_works$x, all_works$y, pch=20, col='#216ea6', cex=0.75)
    points(all_schools$x, all_schools$y, pch=20, col='#dd0000')
    
    #legend('bottomright', legend=c('Workplaces', 'Schools', 'Houses', 'Mosquito paths'),
    #       col=c('#216ea6','#dd0000','#a24f2b','#619e6a80'), pt.cex=c(2, 2.5, 0.5, NA), pch=c(20, 20, 20, NA), lty=c(NA,NA,NA,1), lwd=c(NA,NA,NA,2))
    box()
    if (make_png) dev.off()
    
}

make_yuc_map = function(make_png=F) {
    if (make_png) {
        png('yuc_model-alt-color-PNAS.png', width=2400, height=1600, res=360)
        par(mar=c(3.1, 3.1, 0.1, 1))
    }
    #xlim = c(-90.52197, -87.41819), ylim = c(19.46334,  21.71292)
#    plot(yuc, axes=T, bg='#a6cae0', xlab='Longitude', ylab='Latitude', xlim = c(-90.62197, -87.31819), ylim = c(19.36334,  21.61292))
    plot(yuc, axes=F, bg='#a6cae0', xlab='', ylab='')
    axis(1, mgp=c(3, .5, 0))
    title(xlab='Longitude', line=2)
    axis(2, mgp=c(3, .55, 0))
    title(ylab='Latitude', line=2)
    plot(mex, add=T, col='#999999')
    plot(yuc, add=T, col='#e3debf', border='#00000000')
    plot(mun, add=T, border='white', lwd=0.5)
    plot(yuc, add=T)
    box()

    .xall = c(-90.52197, -87.41819)
    .yall = c(19.46822,  21.70804)
    .xlim = .xall
    .ylim = .yall
    
    pop_movements = get_people_movements(.xlim, .ylim)
    pts = as.data.frame(pop_movements[,c('x1','y1')])
    coordinates(pts) = ~x1+y1
    rast = raster(ncol=ceiling((extent(pts)@xmax - extent(pts)@xmin)/0.00416667),
                  nrow=ceiling((extent(pts)@ymax - extent(pts)@ymin)/0.00416667))
    extent(rast)=extent(pts)
    rast = rasterize(pts, rast, fun='count')
#    pop_colors = colorRampPalette(c("#081d58", "#253494", "#225ea8", "#1d91c0", "#41b6c4", "#7fcdbb", "#c7e9b4", "white"))(100)
    pop_colors = colorRampPalette(c("#253494", "#225ea8", "#1d91c0", "#41b6c4", "#7fcdbb"))(100)
    alpha = 1.0
    pop_colors = paste(pop_colors, sprintf("%x", ceiling(255*alpha)), sep="")
    image(rast, add=T, col=pop_colors)
    #plot(rast, horizontal=F, smallplot=c(.80, .82, .24, .5), legend.only=TRUE, col=pop_colors)

    plot(rast, 
         #horizontal=T, smallplot=c(.08, .45, .93, .95), 
         horizontal=T, smallplot=c(.11, .45, .92, .94), 
         legend.only=TRUE, col=pop_colors,
         #legend.args=list(text='Population density', side=4, font=2, line=2.7))
         legend.args=list(text='Population density', side=3, font=2, line=0.1),
         axis.args=list(mgp=c(3, .3, 0), tcl=-0.3)
    )
    if (make_png) dev.off()
}
