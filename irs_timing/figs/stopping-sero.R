
args <- commandArgs(trailingOnly = TRUE)
# args <- args <- paste0("~/Dropbox/who/fig1_data/stopping-",c("baseline","interventions","sero"),".rds")

require(data.table)

baseline.dt <- readRDS(args[1])
interventions.dt <- readRDS(args[2])

sero.dt <- rbind(
  interventions.dt[,{
      qs <- quantile(seropositive, probs = (0:4)/4)
      q.mn <- mean(seropositive)
      list(q.min = qs[1], q.lo = qs[2], q.med = qs[3], q.hi = qs[4], q.max = qs[5], q.mn=q.mn)
    }, by=.(year, coverage, duration, end_year)
  ],
  baseline.dt[,{
      qs <- quantile(seropositive, probs = (0:4)/4)
      q.mn <- mean(seropositive)
      list(coverage=0, duration=0, end_year=0, q.min = qs[1], q.lo = qs[2], q.med = qs[3], q.hi = qs[4], q.max = qs[5], q.mn=q.mn)
    },by=.(year)
  ]
)

saveRDS(sero.dt, args[3])
