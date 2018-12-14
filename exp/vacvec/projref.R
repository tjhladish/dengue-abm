# reference definitions for project

# in DB, vac_mech = 0/1
vac_names <- c("cmdvi", "traditional", "none")
#' translates vac_mech (and vac = 1 or 0, if needed)
#'
#' @param mech, a 0 (cmdvi) or 1 (traditional)
#' @param vac, a 0 (no vaccination) or 1 (vaccination).  If not supplied, assumed = 1
#' @return ordered factor corresponding to vac_names
trans_vac <- function(mech, vac=rep(1, length(mech))) factor(
  ifelse(vac, vac_names[mech + 1], vac_names[3]),
  levels = vac_names, ordered = T
)

vac_cols <- c("blue", "green", "black")
names(vac_cols) <- vac_names

# in DB, catchup = 0/1
catchup_names <- c("none", "catchup")
#' translates catchup (and vac = 1 or 0, if needed)
#'
#' @param cu, a 0 (none) or 1 (catchup)
#' @param vac, a 0 (no vaccination) or 1 (vaccination).  If not supplied, assumed = 1
#' @return ordered factor corresponding to catchup_names
trans_catchup <- function(cu, vac=rep(1, length(cu))) factor(
  ifelse(vac, catchup_names[cu + 1], catchup_names[1]),
  levels = catchup_names, ordered = T
)

reference.scenario <- data.table(
  vc = 0, vac = 0, vc_coverage = 0,
  vaccine = trans_vac(0, 0), catchup = trans_catchup(0, 0)
)

samplecols <- c("particle", "replicate")