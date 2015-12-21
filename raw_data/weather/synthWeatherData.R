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