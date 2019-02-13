suppressPackageStartupMessages({
  require(data.table)
  require(cowplot)
})

warnnonunique <- function(var, variable, collapse = median) {
  if (length(unique(var)) != 1) warning(sprintf("non unique %s", unique(variable)))
  collapse(var)
}

# debugging args for interactive use
args <- c("figref.rda", "rds/effstats.rds", "fig/fig_2.png")

# expected args:
#  1-3 required: reference_results, interventions_results, effectiveness_stats
#  optional: slice of plot facets
#  required: target plot file
args <- commandArgs(trailingOnly = TRUE)

# load the reference digests
load(args[1])
effstats.dt    <- readRDS(args[2])
tar <- tail(args, 1)

vac.eff <- effstats.dt[variable == "vac.eff", .(
    value = warnnonunique(med, variable),
    vc_coverage = 0,
    scenario = trans_scnario(0, 1)
  ), keyby=.(
    vaccine, catchup = ifelse(catchup=="vc+vac","vac-only","routine"), year
    #, measure = trans_meas(gsub("vac\\.","",variable))
  )
]

vec.eff <- cbind(effstats.dt[variable == "vec.eff", .(
    value = warnnonunique(med, variable),
    scenario = trans_scnario(1, 0)
  ), keyby=.(
    vc_coverage, year
    #, measure = trans_meas(gsub("vec\\.","",variable))
  )
], reference.scenario[,.(vaccine, catchup)])

limits.dt <- rbind(vac.eff, vec.eff)[,
  .(value=c(0, 1), year=-1),
  by=scenario
]

leg.sz <- 0.7

shared <- list(theme_minimal(), aes(
    x=year + 1, y=value, color=scenario,
    fill=interaction(vaccine, catchup), shape=interaction(vaccine, catchup), size=factor(vc_coverage),
    group=interaction(scenario, catchup, vaccine, vc_coverage)
  ), geom_limits(limits.dt),
  scale_color_scenario(guide = "none"),
  scale_year(),
  geom_line(linejoin = "mitre", lineend = "butt"),
  theme(
  	legend.margin = margin(), legend.spacing = unit(25, "pt"),
  	legend.spacing.x = unit(-2,"pt"),
  	legend.key.width = unit(36,"pt"),
  	legend.text = element_text(size=rel(leg.sz)),
  	legend.title = element_text(size=rel(leg.sz)), legend.title.align = 0.5,
  	legend.key.height = unit(1,"pt"),
  	legend.box.spacing = unit(2.5, "pt"),
  	axis.text.y = element_text()
  )
)

pvec <- ggplot(
  vec.eff
) + shared +
scale_size_vectorcontrol(
  breaks=vc_lvls[2:4], guide = gds(1,
    override.aes=list(color=rep(scn_cols["vc"], 3))
  )
)

veclegend <- get_legend(pvec)

#pch.labs <- vac_labels[grep(vacu_lvls)]

pvac <- ggplot(
  vac.eff
) + shared +
	geom_point(size=.75, show.legend = F) +
	geom_point(data=vac.eff[pchstride(year)], size=pchsize) +
  scale_size_vectorcontrol(guide="none") +
	scale_vaccu_interaction(direction = "vertical", label.position="right")

vaclegend <- get_legend(pvac)

basep <- ggplot(
	rbind(vac.eff, vec.eff)
) + facet_grid(scenario ~ ., labeller = facet_labels) + shared +
	theme(
	  axis.title = element_text(size=rel(0.9)),
	  # axis.text = element_text(size=rel(0.7)),
		legend.position = "none",
		panel.spacing.y = unit(15,"pt"),
		strip.text.y = element_text(angle=90),
		plot.margin = margin(t=unit(6,"pt"))
	) +
  geom_pchline(vac.eff) +
	scale_size_vectorcontrol() +
	scale_vaccu_interaction() + 
	scale_effectiveness() + coord_cartesian(clip="off")

p <- ggdraw(basep) + draw_grob(veclegend, x=0.2, y=0.4) + draw_grob(vaclegend, x=0.2, y=-0.07)

save_plot(tar, p, ncol = 1, nrow = 2, base_width = 3.75, base_height = baseh)

# plotutil(p, h=4.5, w=2.75, tar)