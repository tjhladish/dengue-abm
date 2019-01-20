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

cmb.eff <- effstats.dt[variable == "combo.eff", .(
	value = warnnonunique(med, variable),
	scenario = trans_scnario(1, 1)
), keyby=.(vaccine, catchup, vc_coverage, year)]

plot.dt <- rbind(vec.eff, vac.eff)

if (grepl("alt", tar)) plot.dt <- rbind(plot.dt, cmb.eff)

limits.dt <- plot.dt[,
  .(value=c(0, 1), year=-1),
  by=scenario
]

p <- ggplot(
  plot.dt
) + theme_minimal() + aes(
  x=year + 1, y=value, color=scenario,
  fill=catchup, shape=vaccine, linetype=vaccine, size=factor(vc_coverage),
  group=interaction(scenario, catchup, vaccine, vc_coverage)
) +
  facet_grid_freey(scenario ~ ., labeller = facet_labels) +
#  geom_segment(mapping = aes(yend=value, xend=year+1)) +
# geom_step() +
  geom_limits(limits.dt) +
  geom_line(linejoin = "mitre", lineend = "butt") +
  geom_point(size=1) +
  scale_size_vectorcontrol(breaks=vc_lvls[2:4], guide=gds(
    1,
    override.aes=list(
      shape=rep(NA,3 ),
      color=rep(scn_cols["vc"])
    ))
  ) +
	scale_color_scenario(
		name = gsub(" ", "\n", scn_name),
		guide = gds(2, direction="vertical",
			override.aes = list(
				shape=c(vc=vac_pchs["none"], vac=vac_pchs["edv"]),
				fill=scn_cols[c("vc","vac")] # TODO figure out how to make this work?
			)
		)
	) +
  scale_pchlty_vaccine(
  	name = gsub(" ", "\n", vac_name),
  	guide=gds(3, direction="vertical", label.position = "right"), breaks=c("cmdvi", "edv")
  ) +
  scale_fill_catchup(name=gsub(" ", "\n", cu_name),
  	breaks = c("routine", "vac-only"),
  	guide=gds(4, direction="vertical", label.position = "right",
  		override.aes = list(
  			shape = c(21, 21),
  			color = scn_cols[c("vac","vac")]
  		)
  	),
  	na.value=NA
  ) +
  scale_year() +
  scale_y_continuous(expand = c(0,0)) +
  theme(
    legend.margin = margin(), legend.spacing = unit(25, "pt"),
    legend.spacing.x = unit(-2,"pt"),
    legend.text = element_text(size=rel(0.5)),
    legend.title = element_text(size=rel(0.6)), legend.title.align = 0.5,
    panel.spacing.y = unit(15, "pt"), # panel.spacing.x = unit(15, "pt"),
    legend.key.height = unit(1,"pt"),
    legend.box.spacing = unit(2.5, "pt")
  )

plotutil(p, h=7.5, w=3.25, tar)

## TODO SI version:
##  - add cumulative effectiveness row
##  - "explode" interventions into sub levels
##  (two vaccines? maybe 4 vac*catchup & 3 vc coverage)
##  - add IQs
##  - basically whole page, landscape fig