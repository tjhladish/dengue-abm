suppressPackageStartupMessages({
  require(data.table)
  require(cowplot)
  source("utils.R")
  source("labelref.R")
})

warnnonunique <- function(var, variable, collapse = median) {
  if (length(unique(var)) != 1) warning(sprintf("non unique %s", unique(variable)))
  collapse(var)
}

# debugging args for interactive use
args <- c("effstats.rds", "raw_effectiveness_plot.png")

# expected args:
#  1-3 required: reference_results, interventions_results, effectiveness_stats
#  optional: slice of plot facets
#  required: target plot file
args <- commandArgs(trailingOnly = TRUE)

# load the reference digests
effstats.dt    <- readRDS(args[1])

vac.eff <- effstats.dt[variable == "vac.eff", .(
    value = warnnonunique(med, variable),
    vc_coverage=0,
    scenario = trans_scen("vac")
  ), keyby=.(
    vaccine, catchup, year,
    measure = trans_meas(gsub("vac\\.","",variable))
  )
]

vec.eff <- effstats.dt[variable == "vec.eff", .(
    value = warnnonunique(med, variable),
    scenario = trans_scen("vec"),
    vaccine = reference.scenario$vaccine[1],
    catchup = reference.scenario$catchup[1]
  ), keyby=.(
    vc_coverage, year,
    measure = trans_meas(gsub("vec\\.","",variable))
  )
]

combo.eff <- effstats.dt[variable == "combo.eff", .(
  value = warnnonunique(med, variable),
  scenario = trans_scen("comb")
), keyby=.(
  vaccine, catchup, vc_coverage, year,
  measure = trans_meas(gsub("vec\\.","",variable))
)
]

plot.dt <- rbind(vec.eff, vac.eff, combo.eff)

# extend <- copy(plot.dt[year == 39])[, year := 40 ]

gds <- function(ord) guide_legend(
  title.position = "top", direction = "vertical", order = ord,
  label.position = "top"
)

limits.dt <- plot.dt[,
  .(value=c(0, 1), year=-1),
  by=scenario
]

p <- ggplot(
  plot.dt
) + theme_minimal() + aes(
  x=year + 1, y=value,
  color=catchup, linetype=vaccine, size=factor(vc_coverage)
) +
  facet_grid_freey(scenario ~ ., labeller = facet_labels) +
#  geom_segment(mapping = aes(yend=value, xend=year+1)) +
# geom_step() +
  geom_limits(limits.dt) +
  geom_line(linejoin = "mitre", lineend = "butt") +
  scale_linetype_vaccine(guide=gds(2)) +
  scale_color_catchup(guide=gds(3)) +
  scale_size_vectorcontrol(name="Vector Control\nCoverage %", guide=gds(1)) +
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