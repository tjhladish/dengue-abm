
args <- commandArgs(trailingOnly = TRUE)
# args <- paste0("~/Dropbox/who/fig1_data/stopping-",c("baseline","interventions","eff"),".rds")
# args <- c("~/Dropbox/who/fig1_data/foi.baseline.rds", "~/Dropbox/who/fig1_data/foi.interventions.rds", "~/Dropbox/who/fig1_data/foi-eff.rds")

require(data.table)

baseline.dt <- readRDS(args[1])
bkey <- grep("^(imm|s)",names(baseline.dt),invert=T, value=T)
interventions.dt <- readRDS(args[2])
ikey <- grep("^(imm|s|particle)", names(interventions.dt),invert=T, value=T)

join.dt <- interventions.dt[baseline.dt, on=bkey]

eff.dt <- join.dt[!is.na(doy),
  .(eff=ifelse(i.s==s, 0, (i.s-s)/i.s)),
  by=c(ikey, "particle")
]

stat.eff.dt <- eff.dt[,{
    qs <- quantile(eff, probs = (0:4)/4, na.rm = T)
    q.mn <- mean(eff, na.rm = T)
    list(q.min = qs[1], q.lo = qs[2], q.med = qs[3], q.hi = qs[4], q.max = qs[5], q.mn=q.mn)
  },
  by=ikey
]

saveRDS(stat.eff.dt, tail(args,1))
