suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
  source("utils.R")
  source("labelref.R")
})

warnnonunique <- function(var, variable, collapse = median) {
  if (length(unique(var)) != 1) warning(sprintf("non unique %s", unique(variable)))
  collapse(var)
}

# example inputs for debugging
args <- c("baseline.rds","intervention.rds", "effstats.rds", "new_fig_2.png")

# expected args:
#  1-3 required: reference_results, interventions_results, effectiveness_stats
#  optional: slice of plot facets
#  required: target plot file
args <- commandArgs(trailingOnly = TRUE)

# load the reference digests
baseline.dt     <- readRDS(args[1])
intervention.dt <- readRDS(args[2])
effstats.dt     <- readRDS(args[3])

# introspect the scenario keys - all the keys for intervention.dt
#  *except* those specified by `samplecols`
skeys <- setdiff(key(effstats.dt), "variable")

# get the baseline incidence
base.inc <- setkeyv(cbind(baseline.dt[,
  c(dtquantiles(quantile_probs, s),.(scenario = trans_scen("none"))),
  keyby=year
], reference.scenario), skeys)
# cbind ref scenario (no vc, vac) to get right cols for aes

# get the intervention incidence by scenario
inte.inc <- intervention.dt[
  xor(vc, vac), # only want single interventions
  c(dtquantiles(quantile_probs, s), .(scenario = trans_scen(ifelse(vc[1],"vec","vac")))),
  keyby = c(skeys,"vac","vc")
]

# combine ref / scenario incidence,
# scale by yucpopulation, annotate with measure column (for faceting)
inc.dt <- rbind(base.inc, inte.inc)[, 
  .(
    value = med/yucpop,
    measure = trans_meas("inc")
  ), keyby=c("scenario", skeys)
]

# extract the vaccine-only effectiveness measures
# `warnnonunique` should be a no-op with full HPC results
vac.eff <- effstats.dt[variable %in% c("vac.eff","c.vac.eff"), .(
    value = warnnonunique(med, variable),
    vc_coverage = reference.scenario$vc_coverage[1],
    scenario = trans_scen("vac")
  ), keyby=.(
    vaccine, catchup, year,
    measure = trans_meas(gsub("vac\\.","",variable))
  )
]

vec.eff <- effstats.dt[variable %in% c("vec.eff","c.vec.eff"), .(
    value = warnnonunique(med, variable),
    scenario = trans_scen("vec"),
    vaccine = reference.scenario$vaccine[1],
    catchup = reference.scenario$catchup[1]
  ), keyby=.(
    vc_coverage, year,
    measure = trans_meas(gsub("vec\\.","",variable))
  )
]

allplot.dt <- rbind(inc.dt, vec.eff, vac.eff)

extend <- copy(allplot.dt[year == 39])[, year := 40 ]

gds <- function(ord) guide_legend(
  title.position = "top", direction = "horizontal", order = ord,
  label.position = "top"
)

legend_theme <- theme(
  legend.margin = margin(), legend.spacing = unit(25, "pt"),
  legend.text = element_text(size=rel(0.5)),
  legend.title = element_text(size=rel(0.6), margin = margin()),
  legend.title.align = 0.5,
  legend.box.margin = margin(),
  panel.spacing.y = unit(15, "pt"),
  legend.key.height = unit(1,"pt"),
  legend.direction = "horizontal", legend.box = "horizontal",
  legend.justification = c(0.5, 1),
  legend.position = c(0.5, 0),
  axis.title.x = element_text(size=rel(0.7), margin = margin(t=1))
)

limits.dt <- data.table(
  year = -1,
  value=c(
    c(0,1),
    c(0,1),
    c(0, inc.dt[,ceiling(max(value)/300)*300])
  ),
  measure=trans_meas(c("eff","eff","c.eff","c.eff","inc","inc"))
)

ggplot(allplot.dt[
  (scenario == "none") |
  (scenario == "vec" & vc_coverage == 75) |
  catchup == "catchup"
]) + aes(
  x = year+1, y = value,
  color = scenario,
  linetype = vaccine 
) + theme_minimal() + legend_theme +
  facet_grid_freey(measure ~ ., labeller = facet_labels) +
  geom_line() +
  geom_limits(limits.dt) +
  scale_linetype_vaccine(guide=gds(1)) +
  scale_color_scenario(guide=gds(2)) +
  scale_y_continuous(expand = c(0,0)) +
  scale_year()

plotutil(p, h=5, w=7.5, args)