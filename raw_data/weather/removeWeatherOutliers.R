# remove weather outliers
# TODO: should be able to apply same procedure to merida and miami data, and see no removals for miami

require(data.table)
require(reshape2)
merida <- readRDS("merida.RData")
merida.outliers <- melt(merida, id.vars = c("DATE","doy"), variable.name = "extrema", value.name = "celsius")[!is.na(celsius)][,location := factor("MERIDA")][,outlier:=FALSE]

# check weather data visually, tune out outliers until satisfied
merida.outliers[extrema == 'TMAX' & (celsius > 50 | celsius < 15), outlier := TRUE]
merida.outliers[extrema == 'TMIN' & (celsius > 27 | celsius < 5), outlier := TRUE]

merida.trim <- setkey(dcast.data.table(merida.outliers[outlier == FALSE], DATE + doy ~ extrema, value.var = "celsius"), DATE)

saveRDS(merida.trim, "meridaTRIM.RData")

require(ggplot2)

merida.outliers[, year:=year(DATE) ]
pre.extension <- merida.outliers[doy > 330, list(doy = -(-doy %% 365), year = year+1, celsius, outlier, location, extrema)]
post.extension <- merida.outliers[doy < 30, list(doy = doy+365, year = year-1, celsius, outlier, location, extrema)]

extended <- rbind(merida.outliers, pre.extension, post.extension, fill=T)

ggplot(extended) + theme_bw() +
  annotate("rect", xmin=-(-331 %% 365), xmax=0, ymin=-65, ymax=65, fill="grey", alpha=0.1) +
  annotate("rect", xmin=365, xmax=365+29, ymin=-65, ymax=65, fill="grey", alpha=0.1) +
  aes(x=doy, y=celsius, color=extrema, group=interaction(year, extrema)) + geom_line(data=extended[outlier==FALSE], alpha=0.1) + geom_point(mapping=aes(size=outlier)) +
  scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) +
  scale_color_manual(values=c(TMAX='red',TMIN='blue')) + scale_size_manual(breaks=c("TRUE"), values=c(`TRUE`=3,`FALSE`=0))