# make command line output less verbose
suppressPackageStartupMessages({
  require(data.table)
  require(RSQLite)
})

# developer args
args <- c('vacvec-new_yuc.sqlite', "projref.R", "baseline.rds")
args <- c('vacvec-new_yuc.sqlite', "projref.R", "intervention.rds")

# actual args when used with shell
args <- commandArgs(trailingOnly = TRUE)
# args <- c("~/Dropbox/irs_timing-refit0_intro-fix.sqlite", "~/Dropbox/who/fig1_data/baseline.rds")

# db should be first arg
srcdb <- args[1]
# target rds should be last arg
tar <- tail(args, 1)

# TODO convert to rda for shared functions?
source(args[2])
source("utils.R")

# TODO change to input dependency
# samplecols from projref.R: c(particle, replicate)
if (tar == "baseline.rds") {
  qry <- sprintf("select M.*,
          posterior as %s, CAST(realization AS INT) as %s
          from met M, par P, job J
          where M.serial = P.serial and M.serial = J.serial
          and status = 'D'
          and vector_control = 0
          and vac = 0;", samplecols[1], samplecols[2])
} else if (tar == "intervention.rds") {
  qry <- sprintf("select M.*,
          posterior as %s, CAST(realization AS INT) as %s,
          vector_control as vc,
          vac,
          vc_coverage*100 as vc_coverage,
          vac_mech,
          catchup
          from met M, par P, job J
          where M.serial = P.serial and M.serial = J.serial
          and status = 'D'
          and (vector_control = 1 or vac = 1);", samplecols[1], samplecols[2])
} else if (tar == "foi_intervention.rds") {
  qry <- sprintf("select M.*,
          posterior as %s, CAST(realization AS INT) as %s,
          vector_control as vc,
          vac,
          vc_coverage*100 as vc_coverage,
          vac_mech,
          catchup,
          foi
          from met M, par P, job J
          where M.serial = P.serial and M.serial = J.serial
          and status = 'D'
          and (vector_control = 1 or vac = 1);", samplecols[1], samplecols[2])
} else if (tar == "foi_baseline.rds") {
  qry <- sprintf("select M.*,
          posterior as %s, CAST(realization AS INT) as %s, foi
          from met M, par P, job J
          where M.serial = P.serial and M.serial = J.serial
          and status = 'D'
          and vector_control = 0
          and vac = 0;", samplecols[1], samplecols[2])
} else stop(sprintf("don't know the appropriate query for %s",tar))

# dbutil from utils.R
tar.dt <- dbutil(srcdb, qry)

# exclude pre-intervention data & job serial
rmv <- c(grep("s_|imm\\d__",names(tar.dt), value = T), "serial")
# data.table syntax for drop columns
tar.dt <- tar.dt[,.SD,.SDcols=-rmv]

idvs <- samplecols
if (grepl("foi",tar)) idvs <- c("foi", idvs)
if (grepl("intervention", tar)) {
  idvs <- c("vc", "vac", "vc_coverage","vaccine","catchup", idvs)
  tar.dt[vc == 0, vc_coverage := 0]
  # translate vaccine & catchup to standardized lingo for post processing
  tar.dt[, vaccine := trans_vac(vac_mech, vac) ]
  tar.dt[, catchup := trans_catchup(catchup, vac) ]
  tar.dt$vac_mech <- NULL
}

parse.meas.yr <- function(dt) {
  dt[,
    year    := as.integer(gsub("(s|imm\\d_)","", variable))
  ][,
    measure := gsub("(s|imm\\d).+","\\1", variable)
  ]
  dt$variable <- NULL
  return(dt)
}

tar.mlt <- parse.meas.yr(melt.data.table(tar.dt, id.vars = idvs))

tar.dt <- dcast.data.table(tar.mlt,
  as.formula(paste(paste(c(idvs,"year"), collapse=" + "), "measure", sep = " ~ ")),
  value.var = "value"
)

tar.dt[, seropositive := 1-imm0 ]

tar.dt[order(year),
  c.s := cumsum(s),
  by=c(idvs)
]

# from stopping intervention

saveRDS(tar.dt, tar)