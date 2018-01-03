
args <- commandArgs(trailingOnly = TRUE)
# args <- c("~/Dropbox/who/fig1_data/irs_timing-summer_winter-foi.sqlite", "~/Dropbox/who/fig1_data/foi.baseline.rds")
# args <- c("~/Dropbox/who/fig1_data/irs_timing-summer_winter-foi.sqlite", "~/Dropbox/who/fig1_data/foi.interventions.rds")

require(RSQLite)
require(data.table)

isBaseline <- ifelse(grepl("baseline", args[2]), TRUE, FALSE)

keys <- c("particle", "foi")

whereclause <- sprintf('WHERE vector_control == %d AND status == "D"', ifelse(isBaseline,0,1))
keyclause <- paste0(keys[-1], collapse=", ")

if (!isBaseline) {
  keys <- c(keys, "doy")
  keyclause <- paste0(keyclause, ", timing + 1 AS doy")
}

drv = dbDriver("SQLite")
db = dbConnect(drv, args[1], flags=SQLITE_RO)

foi.dt = data.table(dbGetQuery(db, sprintf(
  'SELECT %s, CAST(P.serial / 12 AS INT) AS particle, M.*
   FROM met M
   JOIN par P ON M.serial == P.serial
   JOIN job J ON M.serial == J.serial
   %s;',
   keyclause,
   whereclause
)))

dbDisconnect(db)

rmv <- c(grep("s_|imm\\d__",names(foi.dt), value = T), "serial") # no immunological data for foi?
temp.dt <- foi.dt[,.SD,.SDcols=-rmv]

temp.mlt <- melt.data.table(temp.dt, id.vars = keys)

parse.meas.yr <- function(dt) dt[,
  year    := as.integer(gsub("(s|imm\\d_)","", variable))
][,
  measure := gsub("(s|imm\\d).+","\\1", variable)
][, -"variable", with=F]

res.dt <- dcast.data.table(parse.meas.yr(temp.mlt), ... ~ measure)

# res.dt[, seropositive := 1-imm0 ]

saveRDS(res.dt, args[2])
