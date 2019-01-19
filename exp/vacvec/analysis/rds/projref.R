suppressPackageStartupMessages({
	require(data.table)
})
# reference definitions for project

.args <- commandArgs(trailingOnly = TRUE)

# in DB, vac_mech = 0/1
vac_lvls <- c("cmdvi", "edv", "none")
#' translates vac_mech (and vac = 1 or 0, if needed)
#'
#' @param mech, a 0 (cmdvi) or 1 (edv)
#' @param vac, a 0 (no vaccination) or 1 (vaccination).  If not supplied, assumed = 1
#' @return ordered factor corresponding to vac_lvls
trans_vaccine <- function(mech, vac=rep(1, length(mech))) factor(
  ifelse(vac, vac_lvls[mech + 1], vac_lvls[3]),
  levels = vac_lvls, ordered = T
)
trans_vaccine.data.table <- function(dt) dt[, vaccine := trans_vaccine(vac_mech, vac) ]

scn_lvls <- c("ref", "vc", "vac", "vc+vac")
trans_scnario <- function(vc, vac) factor(
  ifelse(vc & vac, "vc+vac",
  ifelse(vc,       "vc",
  ifelse(vac,      "vac",
                   "ref"
  ))),
  levels = scn_lvls, ordered = T
)
trans_scnario.data.table <- function(dt) dt[, scenario := trans_scnario(vc, vac) ]

int_names <- c("false", "true")
trans_int <- function(combo_gt_assumed) factor(int_names[combo_gt_assumed+1], levels=rev(int_names), ordered = TRUE)

# in DB, catchup = 0/1
cu_lvls <- c("vac-only", "vc+vac", "routine", "none") # c("routine", "catchup", "none")
#' translates catchup (and vac = 1 or 0, if needed)
#'
#' @param cu, a 0 (none) or 1 (catchup)
#' @param vac, a 0 (no vaccination) or 1 (vaccination).  If not supplied, assumed = 1
#' @return ordered factor corresponding to cu_lvls
trans_catchup <- function(catchup, vc, vac) factor(
  ifelse(catchup & vac, # & vac might not be necessary, but guarantees assumption in true branch
    ifelse(vc, "vc+vac", "vac-only"), # has catchup; false branch assumes vac==1
    ifelse(vac, "routine", "none") # no catchup
  ),
  levels = cu_lvls, ordered = T
)
trans_catchup.data.table <- function(dt) dt[, catchup := trans_catchup(catchup, vc, vac)]

reference.scenario <- data.table::data.table(
  vc_coverage = 0, vc = 0, vac = 0,
  vac_mech = 0, catchup = 0
)
trans_catchup.data.table(reference.scenario)
trans_vaccine.data.table(reference.scenario)
trans_scnario.data.table(reference.scenario)

samplecols <- c("particle", "replicate")

quantile_probs <- c(0.025, .25, .5, .75, .975)
names(quantile_probs) <- c("lo.lo","lo","med","hi","hi.hi")

save(list=ls(), file = tail(.args,1))