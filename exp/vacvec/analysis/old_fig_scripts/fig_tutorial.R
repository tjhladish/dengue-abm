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

args <- c("effstats.rds","fig_tutorial.png")
args <- commandArgs(trailingOnly = TRUE)

stat.eff.dt <- readRDS(args[1])
skeys <- key(stat.eff.dt)

stat.eff.dt[, measure := trans_meas(variable) ]

limits <- cbind(data.table(
  measure = trans_meas(c(rep("combo.eff", 2), rep("syn", 2))),
  year = -1,
  med = c(c(0, 1), c(-.15,.25))
), reference.scenario[,.(catchup, vaccine)])

vac_only_reference <- stat.eff.dt[variable == "vac.eff",
  .(measure = trans_meas("combo.eff"),
    med=warnnonunique(med, variable),
    vc_coverage = 0),
  keyby=c(setdiff(skeys, "vc_coverage"))
]

p<-ggplot(stat.eff.dt[measure %in% c("combo.eff")]) +
  facet_grid_freey(measure ~ catchup, labeller=facet_labels) +
  aes(x=year+1, y=med,
    linetype=vaccine, size=factor(vc_coverage), color=catchup) +
  geom_line() +
  geom_line(data = vac_only_reference) +
  geom_limits(limits[measure %in% c("combo.eff")]) +
  scale_linetype_vaccine() +
  scale_size_vectorcontrol() +
  scale_color_catchup(guide="none") +
  scale_year() +
  scale_y_continuous(expand = expand_scale(0, 0)) +
  theme_minimal() + theme(
    legend.direction = "horizontal",
    axis.title.y = element_blank(),
    strip.placement = "outside",
    legend.position = c(0.5, 0.5),
    legend.justification = c(0.5, 0.5),
    legend.margin = margin(), legend.spacing = unit(0, "pt"),
    panel.spacing.y = unit(15, "pt")
  )

plotutil(p, h=5, w=7.5, args)