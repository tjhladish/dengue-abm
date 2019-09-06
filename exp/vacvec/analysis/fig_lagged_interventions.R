suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
  require(cowplot)
})

args <- c("figref.rda", "rds/lag_effstats.rds", "rds/nolag_effstats.rds", "fig/fig_5.png")
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
)[variable %in% c("combo.eff","c.combo.eff") & vaccine == "d70e" & catchup == "vc+vac" & vc_coverage == 75]

ref.combo[, measure := trans_meas(gsub("combo.","", variable, fixed = T)) ][, scenario := "simref" ][, ivn_lag := 0][, obs := "reference" ]

lims <- combo.dt[,.(
  year=-1, med=c(round(min(med)*10)/10, 1)
), by=.(vac_first, measure)]

alphafactor <- 1/2.5

fat <- 4/3

ref.sizes <- c(reference=fat,observed=1/2)*vc_sizes["75"]


pbase <- ggplot() + theme_minimal() + aes(
  x=year+1, y=med, color=scenario, group=factor(ivn_lag),
  linetype=obs, size=obs, alpha=obs
) + adds + facet_grid(vac_first ~ measure, labeller = facet_labels) +
  geom_limits(lims) +
  geom_line(data=ref.combo) +
  geom_line(data=combo.dt[vac_first == 1 & year >= (ivn_lag-1)], mapping=aes(color="vc+vac")) +
  geom_line(data=combo.dt[vac_first == 1 & year < ivn_lag], mapping=aes(color="vac")) +
	geom_pchline(dt=combo.dt[vac_first == 1 & year < ivn_lag], mapping=aes(color="vac"), alpha=1) +
  geom_pchline(dt=combo.dt[vac_first == 1 & year >= ivn_lag], mapping=aes(color="vc+vac"), alpha=1) +
  geom_line(data=combo.dt[vac_first == 0 & year >= (ivn_lag-1)], mapping=aes(color="vc+vac")) +
  geom_line(data=combo.dt[vac_first == 0 & year < ivn_lag], mapping=aes(color="vc")) +
  geom_pchline(dt=combo.dt[vac_first == 0], offset = expression(ivn_lag), mapping=aes(color="vc+vac"), alpha=1) +
  scale_year() + scale_effectiveness(name="Effectiveness", expand = c(0.1, 0, 0, 0)) +
  coord_cartesian(xlim=c(0,20), clip="off") +
  theme(
    legend.direction = "vertical",
    # legend.position = c(0.5,0.5), legend.justification = c(0.5, 0.5),
    legend.text = element_text(size=rel(0.6)),
    legend.title = element_text(size=rel(0.7), vjust = 0),
    legend.title.align = 0.5,
    axis.title = element_text(size=rel(0.9)),
    panel.spacing.y = unit(12,"pt"),
    panel.spacing.x = unit(15,"pt"),
    strip.text = element_text(size=rel(0.8)),
    strip.text.y = element_text(angle=90),
    legend.key.width = unit(20, "pt")
  ) +
  scale_linetype_manual(guide = "none", values = c(reference="21", observed="solid")) +
  scale_size_manual(guide="none", values = ref.sizes) +
  scale_alpha_manual(guide="none", values = c(reference=1, observed=alphafactor))

scn_labels["vc"] <- paste0("75% ", scn_labels["vc"])

pleg <- get_legend(pbase + scale_color_scenario(
  name=gsub(" ","\n", scn_name),
  labels = scn_labels,
  guide=guide_legend(
    override.aes = list(
      shape=c(NA,vac_pchs["d70e"],NA,vac_pchs["d70e"]),
      linetype=c("21","blank","blank","blank"),
      size=c(ref.sizes[1],0,0,0)
    )
)))

pleg2 <- get_legend(pbase + scale_color_scenario(
  name=gsub(" ","\n", scn_name),
  labels = scn_labels,
  guide=guide_legend(
    override.aes = list(
      shape=c(NA,vac_pchs["d70e"],NA,vac_pchs["d70e"]),
      fill=c(NA,scn_cols["vac"],NA,scn_cols["vc+vac"]),
      linetype=c("blank","solid","solid","solid")
    )
)))

p <- ggdraw(pbase + scale_color_scenario(guide="none") + theme(plot.margin = margin(r=72))) +
  draw_grob(pleg, x=0.44, y=0.02) + draw_grob(pleg2, x=0.44, y=0.02)
  
save_plot(tar, p, base_height = baseh*0.55, base_width = 3.75, ncol = 2.5, nrow = 2)