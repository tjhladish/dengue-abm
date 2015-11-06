# fill in data gaps via regression

rm(list = ls())

require(data.table)

miami <- readRDS("miami.RData")
merida <- readRDS("meridaTRIM.RData")

joint <- merida[miami]
joint[is.na(doy), doy := i.doy]
joint$i.doy <- NULL

tmax.mod <- lm(TMAX ~ month(DATE) + TMIN + i.TMAX, data = joint)
tmax.alt <- lm(TMAX ~ month(DATE) + i.TMAX, data = joint)
tmin.mod <- lm(TMIN ~ month(DATE) + TMAX + i.TMIN, data = joint)
tmin.alt <- lm(TMIN ~ month(DATE) + i.TMIN, data = joint)

joint[is.na(TMAX), TMAX.reg := predict(tmax.mod, newdata = .SD)]
joint[is.na(TMIN), TMIN.reg := predict(tmin.mod, newdata = .SD)]
joint[is.na(TMIN) & is.na(TMAX), TMAX.reg := predict(tmax.alt, newdata = .SD)]
joint[is.na(TMIN) & is.na(TMAX), TMIN.reg := predict(tmin.alt, newdata = .SD)]

merida.filled <- joint[,list(doy, TMAX=max(TMAX, TMAX.reg, na.rm=T), TMIN=min(TMIN, TMIN.reg, na.rm=T)),keyby=DATE]

saveRDS(merida.filled, "meridaREG.RData")

require(ggplot2)

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
