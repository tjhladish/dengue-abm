require(data.table)
require(RSQLite)

args <- commandArgs(trailingOnly = TRUE)
# args <- c("~/Dropbox/irs_timing-refit0_intro-fix.sqlite", "~/Dropbox/irs_insecticide-durability-effect_intro-fix.sqlite", "~/Dropbox/who/fig1_data/interventions.rds")

drv = dbDriver("SQLite")

tarcols <- sprintf("s%02d",0:9)

db = dbConnect(drv, args[1], flags=SQLITE_RO)

eff.dt <- data.table(dbGetQuery(db,
  sprintf("SELECT
   timing+1 AS doy, vc_coverage*100 AS coverage,
   campaign_duration AS duration, 90 AS durability,
   posterior AS particle,
   %s
   FROM met M
   JOIN par P ON P.serial = M.serial
   JOIN job J ON J.serial = M.serial
   WHERE status = 'D'
   AND vector_control = 1;", paste0("M.",tarcols,collapse=", ")
  )
))

dbDisconnect(db)

#eff.dt[, particle := floor(serial / 954)]

q <- parse(text=paste0(".(cases10=",paste(tarcols, collapse = "+"),")"))

tmp <- eff.dt[,
  eval(q),
  by=.(doy, coverage, duration, durability, particle)
]
tmp[duration == 2, duration := 365]
tmp[duration == 1, duration := 90]
tmp[duration == 0, duration := 1]

db = dbConnect(drv, args[2], flags=SQLITE_RO)
d2 <- data.table(dbGetQuery(db,
  sprintf("SELECT timing+1 AS doy, vc_coverage*100 AS coverage, 90 AS duration,
  eff_days AS durability, %s
  FROM met M
  JOIN par P ON P.serial = M.serial
  JOIN job J ON J.serial = M.serial
  WHERE status = 'D' and vector_control = 1;", paste0("M.",tarcols,collapse=", "))
))
d2[, particle := 1:.N-1L, by=list(doy, durability)]
d2 <- d2[,
  eval(q),
  by=.(doy, coverage, duration, durability, particle)
]

store <- rbind(tmp, d2)
saveRDS(store, args[3])
