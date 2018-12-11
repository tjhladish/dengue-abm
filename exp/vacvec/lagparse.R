require(data.table)
require(RSQLite)

args <- c('vacvec-new_yuc-ivn_lag.sqlite', "lag_intervention.rds")
args <- commandArgs(trailingOnly = TRUE)
# args <- c("~/Dropbox/irs_timing-refit0_intro-fix.sqlite", "~/Dropbox/who/fig1_data/baseline.rds")

drv = dbDriver("SQLite")
db = dbConnect(drv, args[1], flags=SQLITE_RO)
tar <- args[2]


tar.dt <- data.table(dbGetQuery(db,
   "select M.*,
   posterior as particle, CAST(realization AS INT) as replicate,
   vector_control as vc,
   vac,
   vc_coverage*100 as vc_coverage,
   vac_mech,
   catchup,
   vac_first
   from met M, par P, job J
   where M.serial = P.serial and M.serial = J.serial
   and status = 'D'
   and (vector_control = 1 or vac = 1)"
  ))

dbDisconnect(db)

parse.meas.yr <- function(dt) {
  dt[, 
    year := ifelse(grepl("s_|imm\\d__", variable), -1,1)
  ][,
    measure := gsub("(s|imm\\d).+","\\1", variable)
  ]
  dt[,
    year    := as.integer(gsub("(s|imm\\d_)","", variable))
  ]
  dt$variable <- NULL
  dt
}

# from stopping baseline
#rmv <- c(grep("s_|imm\\d__",names(tar.dt), value = T), "serial")
rmv <- "serial"
tar.dt <- tar.dt[,.SD,.SDcols=-rmv]

idvs <- c("vc", "vac", "vc_coverage","vac_mech","catchup", "vac_first", "particle", "replicate")

tar.mlt <- melt.data.table(tar.dt, id.vars = idvs)
parse.meas.yr(tar.mlt)

tar.dt <- dcast.data.table(tar.mlt,
  as.formula(paste(paste(c(idvs,"year"),collapse=" + "), "measure", sep = " ~ ")),
  value.var = "value"
)

tar.dt[, seropositive := 1-imm0 ]

saveRDS(tar.dt, tar)