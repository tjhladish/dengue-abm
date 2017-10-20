
args <- commandArgs(trailingOnly = TRUE)
# args <- c("~/Dropbox/who/irs_stopping-effect_rerun.sqlite","~/Dropbox/who/fig1_data/stopping-baseline.rds")

# TODO: parse filename a bit to eliminate code duplication?

require(RSQLite)
require(data.table)

drv = dbDriver("SQLite")
db = dbConnect(drv, args[1], flags=SQLITE_RO)

## note, WHERE (i.e., vc_coverage, strat_years, etc) might need to change
## likewise, magic number 12
baseline.dt <- data.table(dbGetQuery(db,
  'select CAST(P.serial / 12 AS INT) AS particle, M.*
   from par P, met M, job J
   where P.serial = M.serial
   and vector_control == 0
   and vc_coverage == 0.75
   and strat_years == 50
   and P.serial = J.serial
   and status = \'D\';'
))

dbDisconnect(db)

rmv <- c(grep("s_|imm\\d__",names(baseline.dt), value = T), "serial")
baseline.dt <- baseline.dt[,.SD,.SDcols=-rmv]

baseline.mlt <- melt.data.table(baseline.dt, id.vars = "particle")

parse.meas.yr <- function(dt) dt[,
  year    := as.integer(gsub("(s|imm\\d_)","", variable))
][,
  measure := gsub("(s|imm\\d).+","\\1", variable)
]

parse.meas.yr(baseline.mlt)

baseline.dt <- dcast.data.table(baseline.mlt, particle + year ~ measure)

baseline.dt[, seropositive := 1-imm0 ]

saveRDS(baseline.dt, args[2])
