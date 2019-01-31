suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

args <- c("figref.rda", "rds/lag_effstats.rds", "rds/effstats.rds", "fig/fig_5.png")
args <- commandArgs(trailingOnly = TRUE)

load(args[1])
stat.eff.dt <- readRDS(args[2])[year < 20]
ref.stat.dt <- readRDS(args[3])[year < 20]
tar <- args[4]

ekeys <- key(stat.eff.dt)

combo.dt <- stat.eff.dt[grepl("combo.eff", variable, fixed = T),
  .(med, obs = "observed"),
  keyby=.(
    vac_first, vc_coverage, vaccine, catchup, year,
    measure = trans_meas(gsub("combo.","", variable, fixed = T))
  )
]

ref.combo <- rbind(
  copy(ref.stat.dt)[, vac_first := 1 ],
  copy(ref.stat.dt)[, vac_first := 0 ]
)[variable %in% c("combo.eff","c.combo.eff") & vaccine == "edv" & catchup == "vc+vac" & vc_coverage == 75]

ref.combo[, measure := trans_meas(gsub("combo.","", variable, fixed = T)) ][, obs := "reference" ]

offsetyr <- 2

p <- ggplot() + theme_minimal() + aes(x=year+1, y=med, color=obs) +
  geom_line(data=ref.combo, size=vc_sizes["75"]) +
  geom_point(data=ref.combo[((year+1) %% 5 == 0) | year == 0], size=2) +
  geom_line(data=combo.dt[vac_first == 1 & year <= offsetyr]) +
  geom_line(data=combo.dt[vac_first == 1 & year >= offsetyr], size=vc_sizes["75"]) +
  geom_point(data=combo.dt[vac_first == 1][((year+1) %% 5 == 0) | year == 0], size=2) +
  
  geom_line(data=combo.dt[vac_first == 0], size=vc_sizes["75"]) +
  geom_point(data=combo.dt[vac_first == 0][((year+1+offsetyr) %% 5 == 0) | year == offsetyr], size=2) +
  
  facet_grid(vac_first ~ measure, labeller = facet_labels) +
  scale_fill_interaction(guide="none") +
  scale_year() +
  scale_effectiveness() +
  coord_cartesian(ylim=c(0.75,1), xlim=c(0,20), clip="off") +
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
    labels=c(reference="Simultaneous Reference", observed="Lagged Result")
  )

# TODO dump shaded area, add intervention annotations, change height aspect

plotutil(p, h=3, w=5.75, tar)