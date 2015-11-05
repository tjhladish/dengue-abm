# remove weather outliers
# TODO: should be able to apply same procedure to merida and miami data, and see no removals for miami
merida <- readRDS("merida.RData")
tmax <- loess(TMAX ~ doy, merida)
tmin <- loess(TMIN ~ doy, merida)
doys <- merida[,list(doy=sort(unique(doy)))]
pred.tmax <- predict(tmax, newdata=doys, se=T)
pred.tmin <- predict(tmin, se=T)

matplot(doys$doy,
  cbind(pred.tmax$fit- qt(0.995,pred.tmax$df)*pred.tmax$se,
        pred.tmax$fit,
        pred.tmax$fit+ qt(0.995,pred.tmax$df)*pred.tmax$se),
  type="l"
)

points(merida$doy, merida$TMAX)

lines(doys$doy, pred.tmax$fit)
lines(unique(sort(merida$doy)), pred.tmax$fit- qt(0.975,pred.tmax$df)*pred.tmax$se, lty=2)
lines(unique(sort(merida$doy)), pred.tmax$fit+ qt(0.975,pred.tmax$df)*pred.tmax$se, lty=2)