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

vec.labxy <- plot.dt[
  year == (25+vc_coverage/25) & scenario == "vc",.(year, value=value+.025),
  by=.(vc_coverage, vaccine, scenario, catchup)
]

vac.labxy <- plot.dt[
  year == 25 & scenario == "vac",.(year, value),
  by=.(vc_coverage, vaccine, scenario, catchup)
]

vac.labxy[
  vaccine == "edv",
  value := value + .05*ifelse(catchup == "routine",1,-1)
]

vac.labxy[
  vaccine == "cmdvi",
  value := value + .05*ifelse(catchup == "routine",-1,1)
]

labxy <- rbind(vec.labxy, vac.labxy)

# labels.dt[labxy, on=.(vc_coverage, vaccine, scenario, catchup)]

if (grepl("alt", tar)) plot.dt <- rbind(plot.dt, cmb.eff)

limits.dt <- plot.dt[,
  .(value=c(0, 1), year=-1),
  by=scenario
]

p <- ggplot(
  plot.dt
) + theme_minimal() + aes(
  x=year + 1, y=value, color=scenario,
  fill=catchup, shape=vaccine, size=factor(vc_coverage),
  group=interaction(scenario, catchup, vaccine, vc_coverage)
) +
  facet_grid_freey(scenario ~ ., labeller = facet_labels) +
#  geom_segment(mapping = aes(yend=value, xend=year+1)) +
# geom_step() +
  geom_limits(limits.dt) +
  geom_line(linejoin = "mitre", lineend = "butt") +
  geom_point(size=1) +
  geom_text(
    aes(x=year+1, y=value, label=label),
    labels.dt[labxy, on=.(vc_coverage, vaccine, scenario, catchup)],
    size = 2
  ) +
  scale_size_vectorcontrol(
    breaks=vc_lvls[2:4], guide="none"
  ) +
	scale_color_scenario(
		guide = "none"
	) +
  scale_shape_vaccine(
  	guide = "none"
  ) +
  scale_fill_catchup(
  	breaks = c("routine", "vac-only"),
  	guide="none",
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