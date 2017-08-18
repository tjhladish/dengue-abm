require(data.table)
require(RSQLite)

args <- commandArgs(trailingOnly = TRUE)
# args <- c("/Users/carlpearson/Dropbox/irs_timing-refit0_intro-fix.sqlite", "~/Dropbox/who/irs_insecticide-durability-effect.sqlite", "~/Dropbox/who/fig1_data/interventions.rds")

drv = dbDriver("SQLite")

db = dbConnect(drv, args[1], flags=SQLITE_RO)

eff.dt <- data.table(dbGetQuery(db,
  "SELECT
   timing+1 AS doy, vc_coverage*100 AS coverage,
   campaign_duration AS duration, 90 AS durability,
   posterior AS particle,
   M.*
   FROM met M
   JOIN par P ON P.serial = M.serial
   JOIN job J ON J.serial = M.serial
   WHERE status = 'D'
   AND vector_control = 1;")
)

dbDisconnect(db)

#eff.dt[, particle := floor(serial / 954)]

tmp <- eff.dt[,
  .(cases10=base::sum(.SD)),
  .SDcols=sprintf("s%02d",0:9),
  by=.(doy, coverage, duration, durability, particle)
]
tmp[duration == 2, duration := 365]
tmp[duration == 1, duration := 90]
tmp[duration == 0, duration := 1]

db = dbConnect(drv, args[2], flags=SQLITE_RO)
d2 <- data.table(dbGetQuery(db,
  'select timing+1 AS doy, vc_coverage*100 AS coverage, 90 AS duration, eff_days AS durability, M.*
  from par P, met M, job J
  where P.serial = M.serial
  and P.serial = J.serial
  and status = \'D\'
  and vector_control = 1;'
))
d2[, particle := 1:.N-1L, by=list(doy, durability)]
d2 <- d2[,
  .(cases10=base::sum(.SD)),
  .SDcols=sprintf("s%02d",0:9),
  by=.(doy, coverage, duration, durability, particle)
]

store <- rbind(tmp, d2)
saveRDS(store, args[3])
