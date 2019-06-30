suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

args <- c("figref.rda", "rds/lag_effstats.rds", "rds/nolag_effstats.rds", "fig/SIfig_5.png")
args <- commandArgs(trailingOnly = TRUE)

load(args[1])
stat.eff.dt <- readRDS(args[2])[year < 20]
ref.stat.dt <- readRDS(args[3])[year < 20]
tar <- args[4]

ekeys <- grep("^(scenario|variable)$", key(stat.eff.dt), value = T, invert = T )

combo.dt <- stat.eff.dt[grepl("combo.eff", variable, fixed = T),
  .(lo, med, hi, obs = "observed", measure = trans_meas(gsub("combo.","", variable, fixed = T))),
  keyby = c(ekeys, "variable")
]

ref.combo <- rbind(
  copy(ref.stat.dt)[, vac_first := 1 ],
  copy(ref.stat.dt)[, vac_first := 0 ]
)[variable %in% c("combo.eff","c.combo.eff") & vaccine == "d70e" & catchup == "vc+vac" & vc_coverage == 75]

ref.combo[, measure := trans_meas(gsub("combo.","", variable, fixed = T)) ][, obs := "reference" ]

lims <- combo.dt[,.(
  year=-1, med=c(floor(min(lo)*10)/10, 1)
), by=.(vac_first, measure)]

ref.combo$ivn_lag <- NULL

fat <- 4/3

p <- ggplot() + theme_minimal() + aes(x=year+1, y=med, color=obs, group=vac_first) +
  geom_line(data=ref.combo, size=vc_sizes["75"]*fat, linetype="21") +
  geom_ribbon(aes(color=NULL, ymax=hi, ymin=lo, fill=obs),
    combo.dt[vac_first == 0 & year < ivn_lag], fill=scn_cols["vc"], alpha=0.5
  ) +
  geom_ribbon(aes(color=NULL, ymax=hi, ymin=lo, fill=obs),
    combo.dt[vac_first == 1 & year < ivn_lag], fill=scn_cols["vac"], alpha=0.5
  ) +
  geom_ribbon(aes(color=NULL, ymax=hi, ymin=lo, fill=obs), combo.dt[year >= (ivn_lag-1)], alpha=0.5) +
  #  geom_line(data=combo.dt[vac_first == 1 & year <= ivn_lag]) +
  geom_line(data=combo.dt[vac_first == 1], size=vc_sizes["75"]/2) +
  geom_line(data=combo.dt[vac_first == 1 & year < ivn_lag], size=vc_sizes["75"]/2, color=scn_cols["vac"]) +
	geom_pchline(dt=combo.dt[vac_first == 1 & year < ivn_lag], color=scn_cols["vac"]) +
	geom_pchline(dt=combo.dt[vac_first == 1 & year >= ivn_lag]) +	
  geom_line(data=combo.dt[vac_first == 0], size=vc_sizes["75"]/2) +
  geom_line(data=combo.dt[vac_first == 0 & year < ivn_lag], size=vc_sizes["75"]/2, color=scn_cols["vc"]) +
	geom_pchline(dt=combo.dt[vac_first == 0], offset = expression(ivn_lag)) +
  geom_limits(lims) +
  ggtitle("Combined Intervention Stagger") +
  facet_grid(measure ~ ivn_lag, labeller = labeller(
    measure = c(eff="Annual Eff.", c.eff="Cumulative Eff."),
    ivn_lag = function(lag) sprintf("%s years", lag)
  )) +
  
#  scale_fill_interaction(guide="none") +
  scale_year() +
  scale_effectiveness(name="Effectiveness", breaks = seq(0,1,by=.25)) +
  coord_cartesian(xlim=c(0,20), clip="off") +
  theme(
    plot.title = element_text(size=rel(0.7), hjust = 0.5, margin=margin()),
    strip.text = element_text(size=rel(0.7))
  ) +
  theme(
    legend.direction = "horizontal",
    legend.position = c(0.5,0.5), legend.justification = c(0.5, 0.5),
    legend.text = element_text(size=rel(0.6)),
    legend.title = element_text(size=rel(0.7), vjust = 0),
    legend.title.align = 0.5,
    panel.spacing.y = unit(20,"pt"),
    panel.spacing.x = unit(15,"pt"),
    strip.text = element_text(size=rel(0.6)),
    strip.text.y = element_text(angle=90)
  ) +
  scale_colour_manual(name=NULL,
    values=c(reference="grey",observed="black"),
    labels=c(reference="Simultaneous Reference", observed="Lagged Result"),
    aesthetics = c("color","fill"), guide="none"
  )

# TODO dump shaded area, add intervention annotations, change height aspect

plotutil(p, h=3, w=5.75, tar)