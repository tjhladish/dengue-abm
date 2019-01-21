suppressPackageStartupMessages({
  require(data.table)
	require(cowplot)
})
# plot normalization functions
# elements re-used when building shared plot elements

.args <- c("utils.R","rds/projref.rda", "figref.rda")
.args <- commandArgs(trailingOnly = TRUE)

source(.args[1])
load(.args[2])

# scales:
#  color = scenario (ref/none, vc-only, vac-only, combination)
#  size = vector control coverage
#  pch = vaccine mechanism (none, dengvaxia, edv)
#  pch fill = catchup (routine vs catchup)

# SCENARIO DIMENSIONING

# scn_lvls <- c("ref", "vc", "vac", "vc+vac") # from projref.R
scn_name <- "Intervention Scenario"
scn_labels <- c("None", "Vector Control Only", "Vaccine Only", "Combination")
scn_cols <- c("grey","blue","darkgreen","darkcyan")
names(scn_labels) <- names(scn_cols) <- names(scn_lvls) <- scn_lvls
scale_color_scenario <- scale_generator(
  "color", scn_name, scn_labels, scn_cols
)

# VC COVERAGE DIMENSIONING

vc_lvls <- seq(from=0,to=75,by=25)
vc_name <- "Vector Control\nCoverage %"
vc_labels <- c("none","25%","50%","75%")
vc_sizes <- seq(from=0.4, by=0.3, length.out = length(vc_labels))
names(vc_labels) <- names(vc_sizes) <- names(vc_lvls) <- vc_lvls
scale_size_vectorcontrol <- scale_generator(
  "size", vc_name, vc_labels, vc_sizes
)

# VACCINE MECH DIMENSIONING

vac_lvls <- c("cmdvi","edv","none")
vac_name <- "Vaccine Model"
vac_labels <- c("Dengvaxia", "D70E", "None")
vac_pchs <- c(24,21,NA) # c(15,16,NA)
vac_ltys <- c("12","61","solid")
names(vac_labels) <- names(vac_pchs) <- names(vac_ltys) <- names(vac_lvls) <- vac_lvls
scale_shape_vaccine <- scale_generator(
  "shape", vac_name, vac_labels, vac_pchs
)
scale_linetype_vaccine <- scale_generator(
	"linetype", vac_name, vac_labels, vac_ltys
)
scale_pchlty_vaccine <- function(...) list(
	scale_shape_vaccine(...), scale_linetype_vaccine(...)
)

# CATCHUP DIMENSIONING

cu_lvls <- c("vac-only", "vc+vac", "routine", "none")
cu_name <- "Vaccine Campaign"
cu_labels <- c("Catchup", "Catchup", "Routine-Only", "None")
cu_fills <- c(scn_cols[c("vac","vc+vac")], "white", NA)
cu_alpha <- c(1, 1, 0.3, 0.3)
names(cu_labels) <- names(cu_fills) <- names(cu_lvls) <- cu_lvls
scale_fill_catchup <- scale_generator(
  "fill", cu_name, cu_labels, cu_fills
)
scale_alpha_catchup <- scale_generator(
  "alpha", cu_name, cu_labels, cu_alpha
)

cuscn_fills <- c(cu_fills, scn_cols)
names(cuscn_fills) <- c(cu_lvls, scn_lvls)


## INTERACTIONS - FOR COMBINATION INTERVENTION PLOTS
# int_names <- c("under", "over")
int_name <- "Interaction"
int_labels <- c("Interfere", "Amplify")
int_fills <- c("red", "blue")
names(int_labels) <- names(int_fills) <- int_names
scale_fill_interaction <- scale_generator(
	"fill", int_name, int_labels, int_fills
)

# 
# ## VACCINE MECHANISM
# 
# vac_title <- "Vaccine"
# vac_cols <- c("blue", "green", "black")
# vac_ltys <- c("dashed", "solid", "dotdash")
# vac_labels <- c("Dengvaxia-like", "EDV", "None")
# names(vac_ltys) <- names(vac_cols) <- names(vac_labels) <- vac_names
# scale_color_vaccine <- scale_generator("color", vac_title, vac_labels, vac_cols)
# scale_linetype_vaccine <- scale_generator("linetype", vac_title, vac_labels, vac_ltys)
# 
# ## CATCHUP
# 
# cu_title <- "Intervention"
# cu_labels <- c("Routine Vac.", "Routine + Catchup", "No Vaccination")
# cu_cols <- c("green", "black", "blue")
# cu_ltys <- c("solid", "dashed", "dotted")
# names(cu_ltys) <- names(cu_labels) <- names(cu_cols) <- catchup_names
# scale_linetype_catchup <- scale_generator("linetype", cu_title, cu_labels, cu_ltys)
# scale_color_catchup <- scale_generator("color", cu_title, cu_labels, cu_cols)
# 
# ## VECTOR CONTROL
# vc_title <- "Vector Control Coverage %"
# sizebase <- 0.5
# sizestep <- 0.25
# vc_sizes <- seq(from=sizebase, by=sizestep, length.out = 4)
# vc_labels <- c("None",25,50,75)
# vc_ltys <- c(3, 2, 5, 1) # c(`25`="dotted",`50`="dashed",`75`="solid")
# names(vc_ltys) <- names(vc_labels) <- names(vc_sizes) <- c(0,25,50,75)
# scale_size_vectorcontrol <- function(
#   name = vc_title, labels = vc_labels,
#   values = vc_sizes, ...
# ) scale_size_manual(
#   name=name, labels=labels,
#   values=values, ...
# )
# scale_linetype_vectorcontrol <- function(
#   name = vc_title, labels = vc_labels,
#   values = vc_ltys, ...
# ) scale_linetype_manual(
#   name=name, labels=labels,
#   values=values, ...
# )
# 
scale_year <- function(...) scale_x_continuous("Year", expand = c(0,0), ...)

yucpop <- 18.17734 # 100ks

meas_names <- c("inc", "eff", "c.eff", "combo.eff", "syn")
trans_meas <- function(meas) factor(meas, levels=meas_names, ordered = T )
meas_labels <- c(
  "Annual Incidence per 100k", "Annual Effectiveness", "Cumulative Effectiveness",
  "Combined Effectiveness", "Interaction"
)
names(meas_labels) <- meas_names


facet_labels <- labeller(
  measure = meas_labels,
  scenario = scn_labels,
  catchup = cu_labels,
  vaccine = vac_labels
)

gds <- function(
  order,
  title.position = "top", direction = "horizontal", label.position = "top",
  ...
) guide_legend(
  title.position = title.position, direction = direction, order = order,
  label.position = label.position, ...
)

labels.dt <- data.table(
  label = c(
    vc_labels[c("25","50","75")],
    paste(
      vac_labels[rep(c("cmdvi","edv"), times=2)],
      cu_labels[rep(c("routine","vac-only"), each=2)],
      sep="\n"
    )
  ),
  scenario = c(
    rep("vc", 3),
    rep("vac", 4)
  ),
  vc_coverage = c(
    c(25, 50, 75),
    rep(vc_lvls[c("0")], 4)
  ),
  vaccine = c(
    rep("none", 3),
    rep(c("cmdvi","edv"), times=2)
  ),
  catchup = c(
    rep("none", 3),
    rep(c("routine","vac-only"), each=2)
  )
)

save(list = ls(), file = tail(.args, 1))