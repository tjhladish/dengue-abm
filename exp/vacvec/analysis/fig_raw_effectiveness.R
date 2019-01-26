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

limits.dt <- plot.dt[,
  .(value=c(0, 1), year=-1),
  by=scenario
]

shared <- list(theme_minimal(), aes(
    x=year + 1, y=value, color=scenario,
    fill=interaction(catchup, vaccine), shape=vaccine, size=factor(vc_coverage),
    group=interaction(scenario, catchup, vaccine, vc_coverage)
  ), geom_limits(limits.dt),
  scale_color_scenario(guide = "none"),
  scale_year(),
  geom_line(linejoin = "mitre", lineend = "butt")
)

pvec <- ggplot(
  vec.eff
) + shared +
  scale_size_vectorcontrol(
    breaks=vc_lvls[2:4], guide = gds(1,
      override.aes=list(color=rep(scn_cols["vc"], 3))
    )
  ) +
  scale_shape_vaccine(
    guide = "none"
  ) +
  scale_fill_catchup(
    guide="none", na.value=NA
  ) +
  scale_y_continuous(
    name="Vector Control-Only Annual Effectiveness",
    expand = c(0,0)
  ) +
  theme(
    legend.margin = margin(), legend.spacing = unit(25, "pt"),
    legend.spacing.x = unit(-2,"pt"),
    legend.position = c(35/40, 0.875),
    legend.justification = c(1,1),
    legend.text = element_text(size=rel(0.5)),
    legend.title = element_text(size=rel(0.6)), legend.title.align = 0.5,
    panel.spacing.y = unit(15, "pt"), # panel.spacing.x = unit(15, "pt"),
    legend.key.height = unit(1,"pt"),
    legend.box.spacing = unit(2.5, "pt")
  )

pvac <- ggplot(
  vac.eff
) + theme_minimal() + shared +
  geom_point(size=1) +
  scale_size_vectorcontrol(
    breaks=vc_lvls[2:4], guide="none"
  ) +
  scale_shape_vaccine(guide=gds(1, direction = "vertical", label.position = "right")) +
  scale_fill_catchup(
    labels = c("routine", "", "vac-only", ""),
    breaks = c("routine", "routine", "vac-only", "vac-only"),
    na.value=NA,
    guide=gds(2, direction = "vertical", label.position = "right")
  ) +
  scale_y_continuous(
    name="Vaccine-Only Annual Effectiveness",
    expand = c(0,0)
  ) +
  theme(
    legend.margin = margin(), #legend.spacing = unit(25, "pt"),
    legend.spacing.x = unit(-2,"pt"),
    legend.text = element_text(size=rel(0.5)),
    legend.title = element_text(size=rel(0.6)), legend.title.align = 0.5,
    panel.spacing.y = unit(15, "pt"), # panel.spacing.x = unit(15, "pt"),
    legend.key.height = unit(1,"pt"),
    legend.box.spacing = unit(2.5, "pt"),
    legend.position = c(35/40, 0.875),
    legend.justification = c(1,1)
  )

plotutil(p, h=7.5, w=3.25, tar)

## TODO SI version:
##  - add cumulative effectiveness row
##  - "explode" interventions into sub levels
##  (two vaccines? maybe 4 vac*catchup & 3 vc coverage)
##  - add IQs
##  - basically whole page, landscape fig