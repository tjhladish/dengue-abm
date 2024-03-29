suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
	require(cowplot)
})
# plot normalization functions
# elements re-used when building shared plot elements

.args <- c("utils.R","rds/projref.rda", "figref.rda")
.args <- commandArgs(trailingOnly = TRUE)

source(.args[1])
load(.args[2])

mainfilter <- expression(foi==1.0 & (vaccine == "d70e" | (false_neg==0.20 & false_pos == 0.05)))

# scales:
#  color = scenario (ref/none, vc-only, vac-only, combination)
#  size = vector control coverage
#  pch = vaccine mechanism (none, dengvaxia, d70e)
#  pch fill = catchup (routine vs catchup)


light_cols <- c("#AAAAFF","#77AA77")
names(light_cols) <- c("vc", "vac")

# SCENARIO DIMENSIONING

# scn_lvls <- c("ref", "vc", "vac", "vc+vac") # from projref.R
names(scn_lvls) <- scn_lvls
scn_name <- "Intervention Scenario"
scn_labels <- c("None", "TIRS Only", "Vaccine Only", "Combination", "Naive Estimate", "Naive Estimate", "Simultaneous\nReference")
scn_cols <- c("yellow","blue","darkgreen","black", light_cols["vc"], light_cols["vac"], "grey")
names(scn_labels) <- names(scn_cols) <- c(scn_lvls, "vc+naive", "vac+naive", "simref")
scale_color_scenario <- scale_generator(
  "color", scn_name, scn_labels, scn_cols
)

# VC COVERAGE DIMENSIONING

vc_lvls <- seq(from=0,to=75,by=25)
vc_name <- "TIRS\nCoverage %"
vc_labels <- c("none","25%","50%","75%")
vc_sizes <- seq(from=0.4, by=0.3, length.out = length(vc_labels))
names(vc_labels) <- names(vc_sizes) <- names(vc_lvls) <- vc_lvls
scale_size_vectorcontrol <- scale_generator(
  "size", vc_name, vc_labels, vc_sizes
)

# VACCINE MECH DIMENSIONING

vac_lvls <- c("cydtdv","d70e","t+cydtdv","none")
vac_name <- "Vaccine Model"
vac_labels <- c("CYD-TDV", "D70E", "CYD-TDV", "None")
vac_pchs <- c(24,21,23,NA) # c(15,16,NA)
vac_nofill_pchs <- c(17,19,18,NA) # c(15,16,NA)
# vac_ltys <- c("12","61","solid")
names(vac_labels) <- names(vac_pchs) <- names(vac_nofill_pchs) <- names(vac_lvls) <- vac_lvls # names(vac_ltys) <-
scale_shape_vaccine <- scale_generator(
  "shape", vac_name, vac_labels, vac_pchs
)
scale_shapenofill_vaccine <- scale_generator(
  "shape", vac_name, vac_labels, vac_nofill_pchs
)
# scale_linetype_vaccine <- scale_generator(
# 	"linetype", vac_name, vac_labels, vac_ltys
# )
# scale_pchlty_vaccine <- function(...) list(
# 	scale_shape_vaccine(...), scale_linetype_vaccine(...)
# )

# CATCHUP DIMENSIONING

cu_lvls <- c("vac-only", "vc+vac", "routine", "none")
cu_name <- "Vaccine Campaign"
cu_labels <- c("Routine w/ Catchup", "Routine w/ Catchup", "Routine-Only", "None")
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

# INTERACTION BETWEEN VACCINE & CATCHUP
vacu_lvls <- with(rbind(expand.grid(head(vac_lvls,-1), cu_lvls[-c(2,4)]),expand.grid(tail(vac_lvls,1), cu_lvls[4])), paste(Var1, Var2, sep="."))
vacu_labels <- with(rbind(expand.grid(head(vac_labels,-1), cu_labels[-c(2,4)]),expand.grid(tail(vac_lvls,1), cu_labels[4])), paste(Var1, Var2, sep=", "))
names(vacu_labels) <- vacu_lvls

scale_vaccu_interaction <- {
	revac <- gsub("^(.+)\\..+$","\\1",vacu_lvls)
	vacu_pchs <- vac_pchs[sort(revac)]
	vacu_pch_labels <- vac_labels[sort(revac)]
	names(vacu_pchs) <- names(vacu_pch_labels) <- vacu_lvls[order(revac)]
	
	recu <- gsub("^.+\\.(.+)$","\\1",vacu_lvls)
	vacu_fills <- cu_fills[sort(recu)]
	vacu_fill_labels <- cu_labels[sort(recu)]
	names(vacu_fills) <- names(vacu_fill_labels) <- vacu_lvls[order(recu)]
	
	scale_pch_vaccu <- scale_generator("shape", vac_name, vacu_pch_labels, vacu_pchs)
	scale_fill_vaccu <- scale_generator("fill", cu_name, vacu_fill_labels, vacu_fills)
	function(...,
		pch.labels=vacu_labels[order(revac)][-5], pch.breaks = names(pch.labels)[c(3,4,5,6,1,2)]
	) list(
		scale_pch_vaccu(breaks=pch.breaks, labels=pch.labels,
			guide=gds(1, override.aes=list(color=scn_cols["vac"], fill=cu_fills[gsub("^.+\\.(.+)$","\\1",pch.breaks[c(1,2,3,4)])], size=pchsize), ...)
		),
		scale_fill_vaccu(breaks=vacu_lvls, guide="none")
	)
}


## INTERACTIONS - FOR COMBINATION INTERVENTION PLOTS
# int_names <- c("under", "over")
int_name <- "Interaction"
int_labels <- c("Interference", "Amplification")
int_fills <- c("red", "blue")
names(int_labels) <- names(int_fills) <- int_names
scale_fill_interaction <- scale_generator(
	"fill", int_name, int_labels, int_fills
)

int_alpha <- 0.25

scale_year <- function(name="Year", ...) scale_x_continuous(
  name=name, expand = c(0,0), ...
)
scale_effectiveness <- function(
  name="Annual Effectiveness",
  breaks=seq(0,1,by=.25), expand = c(0,0), ...
) scale_y_continuous(
  name=name,  breaks=breaks, expand = expand, ...
)

yucpop <- 18.17734 # 100ks

meas_names <- c("inc", "eff", "c.eff", "combo.eff", "syn")
trans_meas <- function(meas) factor(meas, levels=meas_names, ordered = T )
meas_labels <- c(
  "Annual Incidence per 100k", "Annual Effectiveness", "Cumulative Effectiveness",
  "Combined Effectiveness", "Interaction"
)
names(meas_labels) <- meas_names

facet_labels <- labeller(
  vc_coverage = vc_labels,
  measure = meas_labels,
  scenario = scn_labels,
  catchup = cu_labels,
  vaccine = vac_labels,
  vac_first = c(`0`="TIRS First", `1`=paste0(vac_labels["d70e"]," First")),
  foi = function(f) {
    res <- sprintf("%i%%", 100*as.numeric(f))
    res[f=="1"] <- paste0(res[f=="1"], " (Fitted Reference)")
    res
  },
  intervention=c(single="Single Interventions",combined="Combined Interventions")
)

TIRSfacettitle <- list(
  ggtitle("TIRS Coverage"),
  theme(plot.title = element_text(size=rel(1.1), hjust = 0.5, margin=margin()))
)
FOIfacettitle <- list(
  ggtitle("Mosquito Population"),
  theme(
    plot.title = element_text(size=rel(1), hjust = 0.5, margin=margin()),
    strip.text = element_text(size=rel(1))
  )
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
      vac_labels[rep(c("t+cydtdv","d70e"), times=2)],
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
    rep(c("t+cydtdv","d70e"), times=2)
  ),
  catchup = c(
    rep("none", 3),
    rep(c("routine","vac-only"), each=2)
  )
)

pchstride <- function(yr, offset=0, stride=5) (((yr+1) %% stride == 0) & yr > offset) | yr == offset
invpchstride <- function(yr, offset=0, stride=5) !(((yr+1) %% stride == 0) | yr == offset) & yr > offset

pchsize <- 2
smallpch <- 0.5

geom_pchline <- function(dt, offset=0, stride=5, var=expression(year), sz=c(small=smallpch, large=pchsize), show.legend=T, ...) list(
  geom_point(data=dt[invpchstride(eval(var), eval(offset), stride)], size=sz["small"], show.legend = F, ...),
  geom_point(data=dt[pchstride(eval(var), eval(offset), stride)], size=sz["large"], show.legend = show.legend, ...)
)

baseh <- 2.75

save(list = ls(), file = tail(.args, 1))
