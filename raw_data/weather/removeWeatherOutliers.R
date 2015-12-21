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

saveRDS(merida.outliers, "meridapreTRIM.RData")
saveRDS(merida.trim, "meridaTRIM.RData")