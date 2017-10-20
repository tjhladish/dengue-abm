
require(data.table)
require(RSQLite)

args <- commandArgs(trailingOnly = TRUE)
# args <- c("~/Dropbox/irs_timing-refit0_intro-fix.sqlite", "~/Dropbox/who/fig1_data/baseline.rds")

drv = dbDriver("SQLite")
db = dbConnect(drv, args[1], flags=SQLITE_RO)

eff.dt <- data.table(dbGetQuery(db,
  'select M.*, posterior AS particle
  from met M
  join par P on P.serial = M.serial
  join job J on J.serial = M.serial
  where status = \'D\'
  and vector_control = 0;'
))

dbDisconnect(db)
#eff.dt[, particle := floor(serial / 954) ]
store <- eff.dt[,
  .(cases10=base::sum(.SD)),
  .SDcols=sprintf("s%02d",0:9),
  by=particle
]

saveRDS(store, args[2])
