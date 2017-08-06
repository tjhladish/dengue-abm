
args <- commandArgs(trailingOnly = TRUE)
# args <- paste0("~/Dropbox/who/fig1_data/stopping-",c("baseline","interventions","eff"),".rds")

require(data.table)

baseline.dt <- readRDS(args[1])
interventions.dt <- readRDS(args[2])

eff.dt <- interventions.dt[baseline.dt,
  on=c("particle","year")
][,.(eff=ifelse(i.s==s, 0, (i.s-s)/i.s)), by=.(coverage, duration, end_year, particle, year)]

stat.eff.dt <- eff.dt[,{
    qs <- quantile(eff, probs = (0:4)/4, na.rm = T)
    q.mn <- mean(eff, na.rm = T)
    list(q.min = qs[1], q.lo = qs[2], q.med = qs[3], q.hi = qs[4], q.max = qs[5], q.mn=q.mn)
  },
  by=.(coverage, duration, end_year, year)
]

saveRDS(stat.eff.dt, args[3])
