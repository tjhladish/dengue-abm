
require(data.table)
require(lubridate)

args <- commandArgs(trailingOnly = TRUE)
# args <- paste0("~/Dropbox/who/fig1_data/",c("baseline.rds", "interventions.rds", "stat.eff10.rds"))

baseline.dt <- readRDS(args[1])
interventions.dt <- readRDS(args[2])

## perform effectiveness calcs
eff10.dt <- interventions.dt[baseline.dt, on="particle"][,
  # join baseline to interventions on particle basis
  # baseline has *only* particle as key
  .(eff10=(i.cases10-cases10)/i.cases10),
  keyby=.(doy, coverage, duration, durability, efficacy, particle)
]

# take stats across particles; maintains separate results by intervention dimensions
stat.eff10.dt <- eff10.dt[,
  .(med.eff10 = stats::median(eff10)),
  keyby=.(doy, coverage, duration, durability, efficacy)
]

running.mean.smooth.n <- 5
convo.factor <- rep(1, running.mean.smooth.n)/running.mean.smooth.n
# take doy-by-doy statistical results (med.eff10), and smooth with running mean (see ?filter)
stat.eff10.dt[duration != 365, # 365 duration has only one datum, no need (or ability) to smooth
  smooth := as.numeric(filter(
    med.eff10, convo.factor, method = "convolution",
    sides = 2, circular = T
    # want running mean to be centered (sides = 2) and wrap the series (circular = T)
  )),
  by=.(coverage, duration, durability, efficacy)
]
stat.eff10.dt[duration == 365, smooth := med.eff10 ]

saveRDS(stat.eff10.dt,tail(args, 1))
