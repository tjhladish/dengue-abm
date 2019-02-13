suppressPackageStartupMessages({
  require(data.table)
  require(cowplot)
})

args <- c("figref.rda", "rds/lag_effstats.rds", "rds/effstats.rds", "fig/fig_5.png")
args <- commandArgs(trailingOnly = TRUE)

load(args[1])
stat.eff.dt <- readRDS(args[2])[year < 20]
ref.stat.dt <- readRDS(args[3])[year < 20]
tar <- tail(args, 1)

ekeys <- key(stat.eff.dt)

# offsetyr <- 2

#lowylim <- .75

hmult <- 1

if (!grepl("^fig/fig", tar)) {
  lowylim <- .25
  hmult <- 3
  
  vec.dt <- stat.eff.dt[
    variable %in% c("vec.eff", "c.vec.eff") & ((vac_first == 0) | year >= ivn_lag) & med >= lowylim
    ][, measure := trans_meas(gsub("vec.","", variable, fixed = T)) ]
  
  vac.dt <- stat.eff.dt[
    variable %in% c("vac.eff", "c.vac.eff") & ((vac_first == 1) | year >= ivn_lag) & med >= lowylim
    ][, measure := trans_meas(gsub("vac.","", variable, fixed = T)) ]
  
  adds <- list(
    geom_line(data=vec.dt, size=vc_sizes["75"], color=light_cols["vc"]),
    geom_line(data=vac.dt, size=vc_sizes["0"], color=light_cols["vac"]),
    geom_point(data=vac.dt[vac_first == 0][pchstride(year, offset=ivn_lag)], size=pchsize, color=light_cols["vac"]),
    geom_point(data=vac.dt[vac_first == 1][pchstride(year)], size=pchsize, color=light_cols["vac"])
  )
} else adds <- list()

combo.dt <- stat.eff.dt[grepl("combo.eff", variable, fixed = T),
  .(med, obs = "observed"),
  keyby=.(
    ivn_lag, vac_first, vc_coverage, vaccine, catchup, year,
    measure = trans_meas(gsub("combo.","", variable, fixed = T))
  )
]

ref.combo <- rbind(
  copy(ref.stat.dt)[, vac_first := 1 ],
  copy(ref.stat.dt)[, vac_first := 0 ]
)[variable %in% c("combo.eff","c.combo.eff") & vaccine == "edv" & catchup == "vc+vac" & vc_coverage == 75]

ref.combo[, measure := trans_meas(gsub("combo.","", variable, fixed = T)) ][, obs := "reference" ][, ivn_lag := 0]

lims <- combo.dt[,.(
  year=-1, med=c(floor(min(med)*10)/10, 1)
), by=.(vac_first, measure)]

alphafactor <- 1/3.5

p <- ggplot() + theme_minimal() + aes(x=year+1, y=med, color=obs, group=factor(ivn_lag)) +
  adds +
  geom_line(data=ref.combo, size=vc_sizes["75"]*2, linetype = "21") +
  geom_line(data=combo.dt[vac_first == 1], size=vc_sizes["75"]/2, alpha=alphafactor) +
  geom_line(data=combo.dt[vac_first == 1 & year < ivn_lag], size=vc_sizes["75"]/2, color=scn_cols["vac"], alpha=alphafactor) +
	geom_pchline(dt=combo.dt[vac_first == 1 & year < ivn_lag], color=scn_cols["vac"]) +
  geom_pchline(dt=combo.dt[vac_first == 1 & year >= ivn_lag]) +
  geom_line(data=combo.dt[vac_first == 0], size=vc_sizes["75"]/2, alpha=alphafactor) +
  geom_line(data=combo.dt[vac_first == 0 & year < ivn_lag], size=vc_sizes["75"]/2, color=scn_cols["vc"], alpha=alphafactor) +
  geom_pchline(dt=combo.dt[vac_first == 0], offset = expression(ivn_lag)) +
  geom_limits(lims) +
  facet_grid(vac_first ~ measure, labeller = facet_labels) +
  
#  scale_fill_interaction(guide="none") +
  scale_year() +
  scale_effectiveness(name="Effectiveness", expand = c(0.1, 0, 0, 0)) +
  coord_cartesian(xlim=c(0,20), clip="off") +
  theme(
    legend.direction = "vertical",
    # legend.position = c(0.5,0.5), legend.justification = c(0.5, 0.5),
    legend.text = element_text(size=rel(0.6)),
    legend.title = element_text(size=rel(0.7), vjust = 0),
    legend.title.align = 0.5,
    panel.spacing.y = unit(12,"pt"),
    panel.spacing.x = unit(15,"pt"),
    strip.text = element_text(size=rel(0.6)),
    strip.text.y = element_text(angle=90),
    legend.key.width = unit(30, "pt")
  ) +
  scale_colour_manual(name=NULL,
    values=c(reference="grey",observed="black"),
    breaks=c("reference", "observed"),
    labels=c(reference="Simultaneous Reference", observed="Lagged Result"),
    guide=guide_legend(label.position = "top", override.aes = list(
    	linetype = c(reference="21", observed="solid"),
    	shape = c(reference=NA, observed=vac_nofill_pchs["edv"]),
    	size = c(vc_sizes["75"]*2, vc_sizes["75"]/2)
    ))
  )

# TODO dump shaded area, add intervention annotations, change height aspect

save_plot(tar, p, base_height = baseh*0.55, base_width = 3.75, ncol = 2.5, nrow = 2)
#plotutil(p, h=3, w=5.75, tar)