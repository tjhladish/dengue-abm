# check weather data visually

require(ggplot2)
require(data.table)

extremaPlotForm <- function(src) melt(src, id.vars = c("doy","location","DATE"), value.name = "celsius", variable.name = "extrema", na.rm = T)

both <- melt(rbind(
  cbind(readRDS("merida.RData"), location=factor("MERIDA", levels=c("MERIDA","MIAMI"), ordered = T)),
  cbind(readRDS("miami.RData"), location=factor("MIAMI", levels=c("MERIDA","MIAMI"), ordered = T))
), id.vars = c("doy","location","DATE"), value.name = "celsius", variable.name = "extrema", na.rm = T)

ggplot(both) + theme_bw() + facet_grid(. ~ location, scales = "free") +
  aes(x=doy, y=celsius, color=extrema, group=interaction(year(DATE), extrema)) +
  geom_line(alpha=0.1) +
  scale_color_manual(values=c(TMAX='red',TMIN='blue'))

merida.outliers <- readRDS("meridapreTRIM.RData")

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

merida.filled <- readRDS("meridaREG.RData")

molten_weather <- function(src) melt(src, id.vars = c("DATE","doy"), value.name = "celsius", variable.name = "extrema")

ggplot(molten_weather(merida.filled)) + theme_bw() +
  aes(x=doy, y=celsius, color=extrema, group=interaction(year(DATE), extrema)) + geom_line(alpha=0.1) +
  geom_point(data=molten_weather(joint[is.na(TMAX), list(TMAX=TMAX.reg, doy), keyby=DATE]), size=1) +
  geom_point(data=molten_weather(joint[is.na(TMIN), list(TMIN=TMIN.reg, doy), keyby=DATE]), size=1) +
  scale_color_manual(values=c(TMAX='red',TMIN='blue'))

ggplot(molten_weather(merida.filled)) + theme_bw() +
  aes(x=DATE, y=celsius, color=extrema) + geom_line(alpha=0.1) +
  geom_point(data=molten_weather(joint[is.na(TMAX), list(TMAX=TMAX.reg, doy), keyby=DATE]), size=1) +
  geom_point(data=molten_weather(joint[is.na(TMIN), list(TMIN=TMIN.reg, doy), keyby=DATE]), size=1) +
  scale_color_manual(values=c(TMAX='red',TMIN='blue'))

## review climate trends

extremesPlot <- function(src) {
  ggplot(src[,list(maxT=max(celsius), minT=min(celsius)), keyby=list(location, year=year(DATE), extrema)]) +
    theme_bw() + facet_grid(. ~ location, scales = "free") +
    aes(x=year, color=extrema) +
    geom_line(aes(y=maxT, linetype="max")) +
    geom_line(aes(y=minT, linetype="min")) +
    geom_smooth(aes(y=maxT, linetype="max"), method="lm") +
    geom_smooth(aes(y=minT, linetype="min"), method="lm") +
    scale_linetype("annual\nextreme\nextrema") +
    scale_color_manual(values=c(TMAX='red',TMIN='blue'))
}

processedMerida <- extremaPlotForm(readRDS("meridaTRIM.RData")[, location := factor("MERIDA", levels=c("MERIDA","MIAMI"), ordered = T)])
regMerida <- extremaPlotForm(readRDS("meridaREG.RData")[, location := factor("MERIDA", levels=c("MERIDA","MIAMI"), ordered = T)])

extremesPlot(both)
extremesPlot(processedMerida)
extremesPlot(regMerida)

teff.dt <- readRDS("Teff.RData")
ggplot(teff.dt[,list(cumEIR = cumsum(EIRDAY), doy), keyby=year]) +
  theme_bw() +
  aes(x=doy,y=cumEIR,color=year,group=year) +
  geom_line()
