suppressPackageStartupMessages({
  require(data.table)
  require(cowplot)
})

warnnonunique <- function(var, variable, collapse = median) {
  if (length(unique(var)) != 1) warning(sprintf("non unique %s", unique(variable)))
  collapse(var)
}

# debugging args for interactive use
args <- c("labelref.rda", "effstats.rds", "fig/fig_2.png")

# expected args:
#  1-3 required: reference_results, interventions_results, effectiveness_stats
#  optional: slice of plot facets
#  required: target plot file
args <- commandArgs(trailingOnly = TRUE)

# load the reference digests
load(args[1])
effstats.dt    <- readRDS(args[2])

vac.eff <- effstats.dt[variable == "vac.eff", .(
    value = warnnonunique(med, variable),
    vc_coverage=0,
    scenario = trans_scnario(0, 1)
  ), keyby=.(
    vaccine, catchup = ifelse(catchup=="vc+vac","vac-only","routine"), year,
    measure = trans_meas(gsub("vac\\.","",variable))
  )
]

vec.eff <- effstats.dt[variable == "vec.eff", .(
    value = warnnonunique(med, variable),
    scenario = trans_scnario(1, 0),
    vaccine = reference.scenario$vaccine[1],
    catchup = reference.scenario$catchup[1]
  ), keyby=.(
    vc_coverage, year,
    measure = trans_meas(gsub("vec\\.","",variable))
  )
]

plot.dt <- rbind(vec.eff, vac.eff)

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
  geom_point() +
  scale_size_vectorcontrol(breaks=rev(vc_lvls), guide=gds(
    1, direction="vertical",
    override.aes=list(
      shape=c(NA,NA,NA,15),
      color=c("blue","blue","blue","darkgreen")
    ))
  ) +
  scale_color_scenario(guide=gds(2, direction="vertical")) +
  scale_shape_vaccine(guide=gds(3, direction="vertical"), breaks=c("cmdvi", "edv")) +
  scale_fill_catchup(guide=gds(4, direction="vertical"), na.value=NA) +
  scale_year() +
  scale_y_continuous(expand = c(0,0)) +
  theme(
    legend.margin = margin(), legend.spacing = unit(25, "pt"),
    legend.text = element_text(size=rel(0.5)),
    legend.title = element_text(size=rel(0.6)), legend.title.align = 0.5,
    panel.spacing.y = unit(15, "pt"), # panel.spacing.x = unit(15, "pt"),
    legend.key.height = unit(1,"pt")
  )

plotutil(p, h=7.5, w=3.25, args)

## TODO SI version:
##  - add cumulative effectiveness row
##  - "explode" interventions into sub levels
##  (two vaccines? maybe 4 vac*catchup & 3 vc coverage)
##  - add IQs
##  - basically whole page, landscape fig