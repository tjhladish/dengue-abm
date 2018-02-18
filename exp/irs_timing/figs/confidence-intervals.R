
require(data.table)
require(lubridate)
require(RSQLite)

args <- commandArgs(trailingOnly = TRUE)
# args <- c("~/Dropbox/irs_timing-refit0_intro-fix.sqlite", "~/Dropbox/irs_insecticide-durability-effect_intro-fix.sqlite", "~/Dropbox/who/fig1_data/effectiveness.rds")

drv = dbDriver("SQLite")

tarcols <- sprintf("s%02d",0:19)

db = dbConnect(drv, args[1], flags=SQLITE_RO)

baseline.dt <- data.table(dbGetQuery(db, sprintf(
 "select %s, posterior AS particle
  from met M
  join par P on P.serial = M.serial
  join job J on J.serial = M.serial
  where status = 'D'
  and vector_control = 0;",
  paste0("M.",tarcols,collapse=", ")
)))

interventions.dt <- data.table(dbGetQuery(db, sprintf(
 "SELECT
  timing+1 AS doy, vc_coverage*100 AS coverage,
  campaign_duration AS duration, eff_days AS durability,
  posterior AS particle,
  %s
  FROM met M
  JOIN par P ON P.serial = M.serial
  JOIN job J ON J.serial = M.serial
  WHERE status = 'D'
  AND vector_control = 1;",
  paste0("M.",tarcols,collapse=", ")
)))

dbDisconnect(db)

# db = dbConnect(drv, args[2], flags=SQLITE_RO)
# 
# tarcols <- sprintf("s%02d",0:09) # only 10 years of durability scenario currently available
# 
# durability.dt <- data.table(dbGetQuery(db, sprintf(
#  "SELECT timing+1 AS doy, vc_coverage*100 AS coverage, 90 AS duration,
#   eff_days AS durability, %s
#   FROM met M
#   JOIN par P ON P.serial = M.serial
#   JOIN job J ON J.serial = M.serial
#   WHERE status = 'D' and vector_control = 1;",
#   paste0("M.",tarcols,collapse=", ")
# )))
# 
# durability.dt[, particle := 1:.N-1L, by=list(doy, durability)]
# 
# dbDisconnect(db)

baseline.mlt <- melt.data.table(baseline.dt, id.vars = "particle")
baseline.mlt[, year := gsub("s", "", variable) ]
baseline.mlt$variable <- NULL

interventions.mlt <- melt.data.table(
  interventions.dt,
  id.vars = grep("^s", names(interventions.dt), value = T, invert = T)
)
durability.mlt <- melt.data.table(
  durability.dt,
  id.vars = grep("^s", names(durability.dt), value = T, invert = T)
)

interventions.mlt <- rbind(interventions.mlt, durability.mlt)

interventions.mlt[, year := gsub("s", "", variable) ]
interventions.mlt$variable <- NULL

## perform effectiveness calcs
eff.dt <- interventions.mlt[baseline.mlt, on=c("particle","year")][,
  # join baseline to interventions on particle basis
  # baseline has *only* particle as key
  .(numcases=(i.value-value), dencases=i.value),
  keyby=.(doy, coverage, duration, durability, particle, year)
]

eff.dt[,
  cnumcases:=cumsum(numcases),
  by=.(doy, coverage, duration, durability, particle)
]

eff.dt[,
  cdencases:=cumsum(dencases),
  by=.(doy, coverage, duration, durability, particle)
]

result <- eff.dt[,
  .(eff = ifelse(dencases == 0, -Inf, numcases/dencases), ceff = ifelse(cdencases == 0, -Inf, cnumcases/cdencases)),
  keyby=.(doy, coverage, duration, durability, particle, year)
]

stat.result <- result[,{
  pslice <- c(.025, .25, .5, .75, .975)
  subnames <- c("ll","ml","md","mh","hh")
  qeff <- as.list(quantile(eff, probs = pslice, names = FALSE))
  names(qeff) <- paste("eff",subnames,sep=".")
  qceff <- as.list(quantile(ceff, probs = pslice, names = FALSE))
  names(qceff) <- paste("ceff",subnames,sep=".")
  c(qeff, qceff)
}, keyby=.(doy, coverage, duration, durability, year)]

# absslice <- stat.result[coverage == 75 & durability == 90 & doy == 155 & duration == 1]

# take stats across particles; maintains separate results by intervention dimensions

saveRDS(stat.result,args[3])
