require(data.table)
require(knitr)

## get script args; for debugging, uncomment section that follows
args <- commandArgs(trailingOnly = TRUE)
args <- c(
  paste0("~/Dropbox/who/fig1_data/",
         c("foi-baseline.cases","foi-interventions.cases"),
         ".rds"),
  "~/Dropbox/who/fig1_data/tab2.tex"
)

modelpop <- 18.2 # 100k

## read in assorted input data
foibase.dt <- readRDS(args[1])
foiints.dt <- readRDS(args[2])

ref <- foiints.dt[foibase.dt, on=c("particle","foi","year")]

block1 <- ref[between(year, 0, 4),.(
  averted = sum(i.cases - cases),
  bcases = sum(i.cases)
), by=.(foi, doy, particle)]

res1 <- block1[, .(
  averted = as.integer(median(averted)/modelpop),
  ceff = signif(median(ifelse(averted == bcases, 1, averted/bcases)),2)
  ),
  by=.(foi, doy)
][doy==148][,.(averted, ceff, years='1-5'),by=foi]

block2 <- ref[between(year, 5, 9),.(
  averted = sum(i.cases - cases),
  bcases = sum(i.cases)
), by=.(foi, doy, particle)]

res2 <- block2[, .(
  averted = as.integer(median(averted)/modelpop),
  ceff = signif(median(ifelse(averted == bcases, 1, averted/bcases)), 2)
  ), by=.(foi, doy)
][doy==148][,.(averted, ceff, years='6-10'),by=foi]

cat(
  kable(res1[res2, on='foi'], format="latex"),
  file=stdout()
)
