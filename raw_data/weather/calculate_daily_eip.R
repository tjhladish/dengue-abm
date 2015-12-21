rm(list=ls())
require(data.table)

dateRange <- list(start=as.Date("1979/1/1"), end=as.Date("2014/1/1")-1)
merida <- readRDS("meridaREG.RData")

## at this point, slight differences in reconstructed values, due to initial (here) vs later (in orig) pruning by date

# from Chan-Johansson
beta0 = 2.9
betaT = -0.08
vr = 1/4.9

eipfun = function(tempC, beta0, betaT, vr) exp(exp(beta0 + betaT*tempC) + vr/2)

merida[,
  EIPMAX := eipfun(TMAX, beta0, betaT, vr)
][,
  EIPMIN := eipfun(TMIN, beta0, betaT, vr)
]
# WRONG:
# [,
#   EIPAVE := (EIPMIN + EIPMAX)/2
# ][,
#   EIRAVE := 1/EIPAVE
# ]

eirs = merida[,
  list(
    EIRMIN=(1/EIPMIN)/2,
    EIRMAX=(1/EIPMAX)/2
  ),
  keyby=DATE
]

# look backwards to find upper limit on start date
limdate = eirs[dim(eirs)[1]:1,list(datelim=cumsum(EIRMAX+EIRMIN), DATE)][datelim > 1][1,DATE]-1

lookaheads = eirs[between(DATE,dateRange$start,min(limdate, dateRange$end)), {
  cumeir = EIRMAX/2 + EIRMIN # assume initial bite is during the day
  today <- DATE
  slice <- eirs[DATE > today]
  ahead <- 0
  while(cumeir < 1) {
    ahead <- ahead + 1
    day <- slice[ahead]
    add = day$EIRMAX
    cumeir = cumeir + add
    if (cumeir < 1) {
      add = day$EIRMIN
      cumeir = cumeir + add
    }
  }
  del = cumeir - 1
  ahead + (1-del/add)*0.5
}, keyby=DATE]

eipToMu <- function(EIP, vr) log(EIP) - vr/2

eiptoTeff = function(EIP, beta0, betaT, vr) (log(eipToMu(EIP, vr)) - beta0)/betaT

teff.dt <- merida[lookaheads][,
  list(EIP = V1, EIRDAY = (1/EIPMAX + 1/EIPMIN)/2), keyby=list(year=year(DATE), doy)
][,
  list(
    Teff = eiptoTeff(EIP, beta0, betaT, vr),
    mu = eipToMu(EIP, vr),
    EIP, EIRDAY
  ),
  keyby=list(year, doy)
]

#plot(teff.dt[,mean(Teff),keyby=doy])
#points(teff.dt[,eiptoTeff(mean(EIP),beta0,betaT,vr),keyby=doy], col="red")

saveRDS(
  teff.dt,
  "Teff.RData"
)

write.table(teff.dt[, list(mu = eipToMu(mean(EIP),vr)), keyby=doy]$mu, file = "dailyMu.csv", row.names = F, col.names = F)

require(reshape2)
require(ggplot2)

ggplot(melt(teff.dt, id.vars="doy", variable.name="measure")) + theme_bw() +
    facet_grid(measure ~ ., scale="free") +
    aes(x=doy, y=value) +
    geom_line()
