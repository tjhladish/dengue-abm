
args <- commandArgs(trailingOnly = TRUE)
# args <- c("~/Dropbox/who/irs_results-refit2-partial.sqlite","~/Dropbox/who/fig1_data/stopping-interventions.rds")

require(RSQLite)
require(data.table)

drv = dbDriver("SQLite")
db = dbConnect(drv, args[1], flags=SQLITE_RO)

# this provides configuration for which scenario to fetch out strat years 50 for comparison
which_doy <- 1 # which of the doys to use [index-wise] for comparison
short.dt <- data.table(dbGetQuery(db,
  "SELECT DISTINCT
   timing + 1 AS doy
   FROM par P
   JOIN job J ON J.serial = P.serial
   WHERE vector_control = 1
   AND strat_years == 10
   AND status = 'D';"))
keepdoy <- as.integer(short.dt[which_doy])

tarcols <- c(sprintf("s%02d",0:39),sprintf("imm0_%02d",0:39),sprintf("imm1_%02d",0:39),sprintf("imm2_%02d",0:39),sprintf("imm3_%02d",0:39),sprintf("imm4_%02d",0:39))

interventions.dt <- data.table(dbGetQuery(db,sprintf(
  "SELECT vc_coverage*100 AS coverage,
   campaign_duration AS duration,
   timing + 1 AS doy, eff_days AS durability,
   posterior AS particle, strat_years AS end_year,
   %s
   FROM met M
   JOIN par P ON P.serial = M.serial
   JOIN job J ON J.serial = M.serial
   WHERE vector_control = 1
   AND doy = %i
   AND status = 'D';",
   paste0("M.",tarcols,collapse=", "), keepdoy
)))

dbility <- interventions.dt[end_year == 10, unique(durability)]
dration <- interventions.dt[end_year == 10, unique(duration)]
interventions.dt <- interventions.dt[(durability %in% dbility) & (duration %in% dration)]

dbDisconnect(db)

interventions.mlt <- melt.data.table(interventions.dt, id.vars = c("coverage","duration","end_year","durability","particle","doy"))

parse.meas.yr <- function(dt) dt[,
  year    := as.integer(gsub("(s|imm\\d_)","", variable))
][,
  measure := gsub("(s|imm\\d).+","\\1", variable)
]

parse.meas.yr(interventions.mlt)

interventions.dt <- dcast.data.table(interventions.mlt, doy + durability + coverage + duration + end_year + particle + year ~ measure)

interventions.dt[, seropositive := 1-imm0 ]

saveRDS(interventions.dt, tail(args,1))
