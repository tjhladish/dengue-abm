require(data.table)
require(knitr)

args <- commandArgs(trailingOnly = TRUE)
# args <- c(paste0("~/Dropbox/who/fig1_data/stopping-",c("eff","sero"),".rds"),"~/Dropbox/who/fig1_data/tabX.tex")

effectiveness.dt  <- readRDS(args[1])
eff.slice <- effectiveness.dt[end_year == 50]
res.eff <- eff.slice[, .(
  annual.eff = q.med,
  frac.annual = q.med/q.med[1],
  cum.eff = cq.med,
  frac.cum = cq.med/cq.med[1],
  year = year + 1
), by=coverage]

stop.slice <- effectiveness.dt[end_year == 10]
stop.eff <- stop.slice[, .(
  annual.eff = q.med,
  frac.annual = q.med/q.med[1],
  cum.eff = cq.med,
  frac.cum = cq.med/cq.med[1],
  year = year + 1
), by=coverage]

seroprevalence.dt <- readRDS(args[2])
sero.slice <- seroprevalence.dt[end_year == 50]
res.sero <- sero.slice[, .(
  sero = q.med,
  year = year + 1
), by=coverage]

res <- res.eff[res.sero, on=c("year","coverage")]

res[,year[which.min(sero)],by=coverage]

cat(
  kable(res.eff[res.sero, on=c("year","coverage")][,.(
    Year=year,
    `Effectiveness`=round(annual.eff,2), `% Initial`=round(frac.annual,2),
    `Effectiveness`=round(cum.eff,2), `% Initial`=round(frac.cum,2),
    `Seropositivity`=round(sero,2)
  )], format="latex"),
  file=stdout()
)
