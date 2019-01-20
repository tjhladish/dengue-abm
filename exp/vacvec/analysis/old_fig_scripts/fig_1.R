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

# example inputs for debugging
args <- c("baseline.rds","intervention.rds","effstats.rds", "only=eff", "fig_1.png")
args <- c("baseline.rds","intervention.rds","effstats.rds", "fig_1.png")

# expected args:
#  1-3 required: reference_results, interventions_results, effectiveness_stats
#  optional: slice of plot facets
#  required: target plot file
args <- commandArgs(trailingOnly = TRUE)

# load the reference digests
baseline.dt     <- readRDS(args[1])
intervention.dt <- readRDS(args[2])
effstats.dt    <- readRDS(args[3])

# introspect the scenario keys - all the keys for intervention.dt
#  *except* those specified by `samplecols` (since we will take stats over those in next
# two steps)
skeys <- setdiff(key(intervention.dt), samplecols)

# get the quantiles for reference incidence
# need to bind the reference scenario (no vc, no vac) to
# facet correctly
base.inc <- setkeyv(cbind(baseline.dt[,
  dtquantiles(quantile_probs, s),
  by=year
], reference.scenario), skeys)

# same, but for interventions.  no need to add scenario info,
# as that is included in skeys
inte.inc <- intervention.dt[
  xor(vc, vac), # consider only single interventions
  dtquantiles(quantile_probs, s),
  keyby = skeys
]

inc.dt <- rbind(
  copy(base.inc)[, scenario := trans_scen("vac") ],
  copy(base.inc)[, scenario := trans_scen("vec") ],
  inte.inc[, scenario := trans_scen(ifelse(vc,"vec","vac")) ]
)[,
  .(value = med/yucpop, measure = trans_meas("inc")),
  keyby = c("scenario", setdiff(skeys,c("vac","vc")))
]

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

twoplot <- rbind(inc.dt, vec.eff, vac.eff)

extend <- copy(twoplot[year == 39])[, year := 40 ]

gds <- function(ord) guide_legend(
  title.position = "top", direction = "horizontal", order = ord,
  label.position = "top"
)

limits.dt <- twoplot[,.(value=c(0, ceiling(max(value, na.rm = T))), year=-1), by=measure]

limits.dt[measure == "inc", value := { ceiling(value/300)*300 }]

if (length(only <- grep("only=", args, value = T))) {
  if (!grepl("only=(eff|inc|vac|vec)", only)) stop(sprintf("Provided incorrect 'only=' arg: %s", only))
  only <- gsub("only=", "", only)
  if ((only == "eff") | (only == "inc")) {
    limits.dt <- limits.dt[grepl(only, measure)]
    filter.dt <- twoplot[grepl(only, measure)]
    legend_theme <- theme(
      legend.position = "bottom"#, legend.justification = c(0.5, 0.5)
    )
  } else {
    filter.dt <- twoplot[grepl(only, scenario)]
    legend_theme <- theme(
      legend.position = "right",
      legend.justification = c(0, 0.5),
      legend.box = "vertical"
    )
  }
} else {
  filter.dt <- twoplot
  legend_theme <- theme(
    legend.position = c(0.5, 0.5),
    legend.justification = c(0.5, 0.5),
    legend.box = "horizontal"
  )
}

p <- ggplot(
#  rbind(twoplot, extend) # for step / segment versions
  filter.dt
) + theme_minimal() + aes(
  x=year + 1,
  y=value, color=catchup, linetype=vaccine,
    size=factor(vc_coverage)
  ) +
  facet_grid_freey(measure ~ scenario, labeller = facet_labels) +
#  geom_segment(mapping = aes(yend=value, xend=year+1)) +
# geom_step() +
  geom_limits(limits.dt) +
  geom_line(linejoin = "mitre", lineend = "butt") +
  scale_linetype_vaccine(guide=gds(1)) +
  scale_color_catchup(guide=gds(2)) +
  scale_size_vectorcontrol(guide=gds(3)) +
  scale_year() +
  scale_y_continuous(expand = c(0,0)) +
  legend_theme +
  theme(
    legend.margin = margin(), legend.spacing = unit(25, "pt"),
    legend.text = element_text(size=rel(0.5)),
    legend.title = element_text(size=rel(0.6)), legend.title.align = 0.5,
    panel.spacing.y = unit(30, "pt"), panel.spacing.x = unit(15, "pt"),
    legend.key.height = unit(1,"pt")
  )

plotutil(p, h=5, w=7.5, args)

## TODO SI version:
##  - add cumulative effectiveness row
##  - "explode" interventions into sub levels
##  (two vaccines? maybe 4 vac*catchup & 3 vc coverage)
##  - add IQs
##  - basically whole page, landscape fig