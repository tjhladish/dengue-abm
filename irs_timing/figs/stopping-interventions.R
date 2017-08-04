
args <- commandArgs(trailingOnly = TRUE)
# args <- c("~/Dropbox/who/irs_stopping-effect_rerun.sqlite","~/Dropbox/who/fig1_data/stopping-interventions.rds")

require(RSQLite)
require(data.table)

drv = dbDriver("SQLite")
db = dbConnect(drv, args[1], flags=SQLITE_RO)

# may need to tweak magic number 12
interventions.dt <- data.table(dbGetQuery(db,
  "SELECT vc_coverage*100 AS coverage,
   campaign_duration AS duration,
   strat_years AS end_year,
   CAST(P.serial / 12 AS INT) AS particle,
   M.*
   FROM par P, met M, job J
   WHERE P.serial = M.serial
   AND vector_control = 1
   AND P.serial = J.serial
   AND status = 'D';"
))

dbDisconnect(db)

rmv <- c(grep("s_|imm\\d__",names(interventions.dt), value = T), "serial")
interventions.dt <- interventions.dt[,.SD,.SDcols=-rmv]

interventions.mlt <- melt.data.table(interventions.dt, id.vars = c("coverage","duration","end_year","particle"))

parse.meas.yr <- function(dt) dt[,
  year    := as.integer(gsub("(s|imm\\d_)","", variable))
][,
  measure := gsub("(s|imm\\d).+","\\1", variable)
]

parse.meas.yr(interventions.mlt)

interventions.dt <- dcast.data.table(interventions.mlt, coverage + duration + end_year + particle + year ~ measure)

interventions.dt[, seropositive := imm1+imm2+imm3+imm4 ]

saveRDS(interventions.dt, args[2])
