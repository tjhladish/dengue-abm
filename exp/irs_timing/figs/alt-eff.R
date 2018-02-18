
args <- commandArgs(trailingOnly = TRUE)
# args <- paste0("~/Dropbox/who/fig1_data/stopping-",c("baseline","interventions","eff"),".rds")
# args <- paste0("~/Dropbox/who/fig1_data/foi-",c("baseline","interventions","eff"),".rds")

require(data.table)

baseline.dt <- readRDS(args[1])
bkey <- grep("^(imm|s)",names(baseline.dt),invert=T, value=T)
interventions.dt <- readRDS(args[2])
ikey <- grep("^(imm|s|particle|year)", names(interventions.dt),invert=T, value=T)

join.dt <- interventions.dt[
  baseline.dt, on=bkey ][,
  .(intcases = s, basecases = i.s), keyby=c(ikey, "particle", "year")
]

calc.dt <- join.dt[,
  .(averted = basecases-intcases, denom = basecases, caverted = cumsum(basecases-intcases), cdenom=cumsum(basecases), year),
  by=c(ikey, "particle")
]

eff.dt <- calc.dt[!is.na(doy),
  .(
    eff=ifelse(averted==0, 0, averted/denom),
    ceff=ifelse(caverted==0, 0, caverted/cdenom),
    year
  ),
  by=c(ikey, "particle")
]

stat.eff.dt <- eff.dt[,{
    qs <- quantile(eff, probs = (0:4)/4, na.rm = T)
    q.mn <- mean(eff, na.rm = T)
    cqs <- quantile(ceff, probs = (0:4)/4, na.rm = T)
    cq.mn <- mean(ceff, na.rm = T)
    list(
      q.min = qs[1], q.lo = qs[2], q.med = qs[3], q.hi = qs[4], q.max = qs[5], q.mn=q.mn,
      cq.min = cqs[1], cq.lo = cqs[2], cq.med = cqs[3], cq.hi = cqs[4], cq.max = cqs[5], cq.mn=cq.mn
    )
  },
  by=c(ikey,"year")
]

saveRDS(stat.eff.dt, tail(args,1))
