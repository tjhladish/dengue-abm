
require(data.table)
require(RSQLite)

args <- commandArgs(trailingOnly = TRUE)

drv = dbDriver("SQLite")
db = dbConnect(drv, args[1], flags=SQLITE_RO)

eff.dt <- data.table(dbGetQuery(db,
  'select M.*
  from par P, met M, job J
  where P.serial = M.serial
  and P.serial = J.serial
  and status = \'D\'
  and vector_control = 0;'
))

dbDisconnect(db)
eff.dt[, particle := floor(serial / 954) ]
store <- eff.dt[,
  .(cases10=base::sum(.SD)),
  .SDcols=sprintf("s%02d",0:9),
  by=particle
]

saveRDS(store, args[2])
