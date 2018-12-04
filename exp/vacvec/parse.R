require(data.table)
require(RSQLite)

args <- c('~/Dropbox/who/vacvec/vacvec.sqlite', "baseline.rds")
args <- c('vacvec-new_yuc.sqlite', "intervention.rds")
args <- commandArgs(trailingOnly = TRUE)
# args <- c("~/Dropbox/irs_timing-refit0_intro-fix.sqlite", "~/Dropbox/who/fig1_data/baseline.rds")


drv = dbDriver("SQLite")
db = dbConnect(drv, args[1], flags=SQLITE_RO)
tar <- args[2]

if (tar == "baseline.rds") {
  tar.dt <- data.table(dbGetQuery(db,
     "select M.*,
     posterior as particle, CAST(realization AS INT) as replicate
     from met M, par P, job J
     where M.serial = P.serial and M.serial = J.serial
     and status = 'D'
     and vector_control = 0
     and vac = 0;"
  ))
} else {
  tar.dt <- data.table(dbGetQuery(db,
     "select M.*,
     posterior as particle, CAST(realization AS INT) as replicate,
     vector_control as vc,
     vac,
     vc_coverage*100 as vc_coverage,
     vac_mech,
     catchup
     from met M, par P, job J
     where M.serial = P.serial and M.serial = J.serial
     and status = 'D'
     and (vector_control = 1 or vac = 1)"
  ))
}

dbDisconnect(db)

parse.meas.yr <- function(dt) {
    dt[,
        year    := as.integer(gsub("(s|imm\\d_)","", variable))
    ][,
        measure := gsub("(s|imm\\d).+","\\1", variable)
    ]
  dt$variable <- NULL
  dt
}

# from stopping baseline
rmv <- c(grep("s_|imm\\d__",names(tar.dt), value = T), "serial")
tar.dt <- tar.dt[,.SD,.SDcols=-rmv]

idvs <- c("particle", "replicate")
if (tar != "baseline.rds") {
  idvs <- c("vc", "vac", "vc_coverage","vac_mech","catchup", idvs)
}

tar.mlt <- melt.data.table(tar.dt, id.vars = idvs)
parse.meas.yr(tar.mlt)

tar.dt <- dcast.data.table(tar.mlt,
  as.formula(paste(paste(c(idvs,"year"),collapse=" + "), "measure", sep = " ~ ")),
  value.var = "value"
)

tar.dt[, seropositive := 1-imm0 ]

# from stopping intervention

saveRDS(tar.dt, tar)