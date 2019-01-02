# plot normalization functions
# elements re-used when building shared plot elements

source("projref.R")

scale_generator <- function(scale, defname, deflabels, defvalues) {
  scalef <- switch(scale,
    color = scale_color_manual,
    fill = scale_fill_manual,
    linetype = scale_linetype_manual,
    size = scale_size_manual
  )
  return(function(
      name=defname, labels=deflabels, values=defvalues,
      ...
  ) scalef(
    name = name, labels = labels,
    values = values, ...
  ))
}

## INTERACTIONS - FOR COMBINATION INTERVENTION PLOTS

int_title <- "Interaction"
int_labels <- c("Interference", "Synergy")
int_cols <- c("red", "blue")
names(int_labels) <- names(int_cols) <- int_names

## VACCINE MECHANISM

vac_title <- "Vaccine"
vac_cols <- c("blue", "green", "black")
vac_ltys <- c("dashed", "solid", "dotdash")
vac_labels <- c("Dengvaxia-like", "EDV", "None")
names(vac_ltys) <- names(vac_cols) <- names(vac_labels) <- vac_names
scale_color_vaccine <- scale_generator("color", vac_title, vac_labels, vac_cols)
scale_linetype_vaccine <- scale_generator("linetype", vac_title, vac_labels, vac_ltys)

## CATCHUP

cu_title <- "Intervention"
cu_labels <- c("Routine Vac.", "Routine + Catchup", "No Vaccination")
cu_cols <- c("green", "black", "blue")
cu_ltys <- c("solid", "dashed", "dotted")
names(cu_ltys) <- names(cu_labels) <- names(cu_cols) <- catchup_names
scale_linetype_catchup <- scale_generator("linetype", cu_title, cu_labels, cu_ltys)
scale_color_catchup <- scale_generator("color", cu_title, cu_labels, cu_cols)

## VECTOR CONTROL
vc_title <- "Vector Control Coverage %"
sizebase <- 0.5
sizestep <- 0.25
vc_sizes <- seq(from=sizebase, by=sizestep, length.out = 4)
vc_labels <- c("None",25,50,75)
vc_ltys <- c(3, 2, 5, 1) # c(`25`="dotted",`50`="dashed",`75`="solid")
names(vc_ltys) <- names(vc_labels) <- names(vc_sizes) <- c(0,25,50,75)
scale_size_vectorcontrol <- function(
  name = vc_title, labels = vc_labels,
  values = vc_sizes, ...
) scale_size_manual(
  name=name, labels=labels,
  values=values, ...
)
scale_linetype_vectorcontrol <- function(
  name = vc_title, labels = vc_labels,
  values = vc_ltys, ...
) scale_linetype_manual(
  name=name, labels=labels,
  values=values, ...
)

scale_year <- function(...) scale_x_continuous("Year", expand = c(0,0), ...)

yucpop <- 18.17734 # 100ks

meas_names <- c("inc", "eff", "c.eff", "combo.eff", "syn")
trans_meas <- function(meas) factor(meas, levels=meas_names, ordered = T )
meas_labels <- c(
  "Annual Incidence per 100k", "Annual Effectiveness", "Cumulative Effectiveness",
  "Combined Effectiveness", "Interaction"
)
names(meas_labels) <- meas_names

scen_title <- "Scenario"
scen_names <- c("none","vec","vac","comb")
trans_scen <- function(scn) factor(scn, levels = scen_names, ordered = T)
scen_labels <- c("Reference", "Vector Control Only", "Vaccine Only", "Combined Interventions")
scen_cols <- c("black", "green", "blue", "purple")
names(scen_cols) <- names(scen_labels) <- scen_names
scale_color_scenario <- function(
  name = scen_title, labels = scen_labels,
  values = scen_cols, ...
) scale_color_manual(
  name = name, labels = labels,
  values = values, ...
)

facet_labels <- labeller(
  measure = meas_labels,
  scenario = scen_labels,
  catchup = cu_labels
)

gds <- function(
  order,
  title.position = "top", direction = "horizontal",
  ...
) guide_legend(
  title.position = title.position, direction = direction, order = order,
  label.position = "top", ...
)