require(data.table)
require(RSQLite)

args <- commandArgs(trailingOnly = T)
# args <- c("~/Dropbox/who/fig1_data/irs_timing-summer_winter-foi.sqlite", "~/Dropbox/who/fig1_data/foi_data.rds")

tarcols <- sprintf("s%02d",0:9)
q <- parse(text=paste0(".(cases10=",paste(tarcols, collapse = "+"),")"))

drv = dbDriver("SQLite")
db  = dbConnect(drv, args[1], flags=SQLITE_RO)

foi.baseline = data.table(dbGetQuery(db, sprintf(
  'SELECT foi, CAST(M.serial/12 as INT) AS particle, %s
   FROM met M
   JOIN par P ON M.serial == P.serial
   JOIN job J ON M.serial == J.serial
   WHERE vector_control == 0
   AND status == "D";', paste0("M.", tarcols, collapse=", ")
)))

b.cases <- foi.baseline[,
  q,
#  .(cases10=base::sum(.SD)),
#  .SDcols=sprintf("s%02d",0:9),
  keyby=.(foi, particle)
]

# foi.baseline[, .N, by=foi][,all(N==1000)]

foi.intervention = data.table(dbGetQuery(db,
  'SELECT timing, foi, CAST(M.serial/12 as INT) AS particle, %s
   FROM met M
   JOIN par P ON M.serial == P.serial
   JOIN job J ON M.serial == J.serial
   WHERE vector_control == 1
   AND status == "D";', paste0("M.", tarcols, collapse=", ")
))

i.cases <- foi.intervention[,
  q,
#  .(cases10=base::sum(.SD)),
#  .SDcols=sprintf("s%02d",0:9),
  keyby=.(foi, particle, timing)
]

dbDisconnect(db)

res <- i.cases[b.cases][,
  .(eff = ifelse(cases10==i.cases10,0,(i.cases10-cases10)/i.cases10)),
  by=.(foi, particle, timing)
][,
  .(med.eff=median(eff)),
  keyby=.(foi,timing)
]

saveRDS(res, args[2])
