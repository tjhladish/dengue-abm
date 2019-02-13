suppressPackageStartupMessages({
  require(data.table)
  require(cowplot)
})

args <- c("figref.rda", "rds/foi_effstats.rds", "rds/effstats.rds", "fig/fig_5.png")
args <- commandArgs(trailingOnly = TRUE)

tar <- tail(args, 1)

load(args[1])

base.stat.eff.dt <- readRDS(args[3])[vaccine == "cmdvi" & catchup == "routine" & vc_coverage == 75]
base.stat.eff.dt[, foi := 1.0 ]

stat.eff.dt <- rbind(readRDS(args[2]), base.stat.eff.dt)[variable %in% c("combo.eff","vac.eff","vec.eff","ind.eff")]
stat.eff.dt[,
  scenario := factor(ifelse(variable == "combo.eff", "vc+vac",
              ifelse(variable == "vac.eff"  , "vac",
              ifelse(variable == "vec.eff"  , "vc", "not"))),  levels = c("vc","vac","vc+vac","not"), ordered = T)
]
stat.eff.dt[scenario == "vac", vc_coverage := 0]
stat.eff.dt[scenario == "vc", vaccine := "none"]

real.dt <- stat.eff.dt[variable %in% c("combo.eff","vac.eff","vec.eff")]

leg.sz <- 0.8

pbase <- ggplot(
  real.dt
) + theme_minimal() +
  aes(shape=vaccine, color=scenario, x=year+1, y=med, size=factor(vc_coverage), fill=catchup) +

  geom_line(mapping=aes(color="vc+naive"), data=stat.eff.dt[variable == "ind.eff"]) +
  geom_pchline(dt=stat.eff.dt[variable == "ind.eff"], color=light_cols["vac"]) +
  geom_line(data=real.dt[scenario == "vc"]) +

  geom_line(data=real.dt[scenario != "vc"]) +
  geom_pchline(dt=real.dt, fill="white") +
  scale_size_vectorcontrol(guide="none") +
  scale_shape_vaccine(guide="none") +
  scale_year() +
  scale_effectiveness() +
	scale_fill_catchup(guide="none") +
  facet_grid(. ~ foi, labeller=facet_labels) +
  FOIfacettitle +
  coord_cartesian(clip="off", ylim=c(-.125,1), xlim=c(0,40)) + theme(
    panel.spacing.x = unit(12, "pt")
  ) + theme(
  	axis.title = element_text(size=rel(1)),
  	axis.text = element_text(size=rel(0.9)),
  	legend.margin = margin(), legend.spacing = unit(25, "pt"),
  	legend.text = element_text(size=rel(leg.sz)),
  	legend.title = element_text(size=rel(leg.sz)), legend.title.align = 0.5,
  	panel.spacing.x = unit(15, "pt"),
  	strip.text = element_text(size=rel(1)),
  	strip.text.y = element_text(angle=90),
  	legend.key.height = unit(1,"pt"),
  	legend.box.spacing = unit(2.5, "pt")
  )

labs <- scn_labels
labs["vc"] <- paste0(vc_labels["75"]," ",scn_labels["vc"])
labs["vac"] <- paste0("Routine ", vac_labels["cmdvi"]," Only")
labs["vc+vac"] <- paste0("TIRS"," & ",vac_labels["cmdvi"])
labs["vc+naive"] <- labs["vac+naive"] <- paste0("Naive ",labs["vc+vac"])

ppchleg <- get_legend(pbase + scale_color_scenario(labels=labs, guide=guide_legend(
	override.aes = list(
		shape = c(vac_pchs["cmdvi"],NA,vac_pchs["cmdvi"],vac_pchs["cmdvi"]),
		linetype = 0,
		color = c(scn_cols["vac"], NA, scn_cols["vac+naive"], scn_cols["vc+vac"])
	)
)))

pltyleg <- get_legend(pbase + scale_color_scenario(labels=labs, guide=guide_legend(
	override.aes = list(
		shape = NA,
		linetype = "solid",
		size = c(vc_sizes["0"], rep(vc_sizes["75"], 3))
	)
)))

p <- ggdraw(pbase + scale_color_scenario(guide="none")) + draw_grob(pltyleg, x=0.05, y=.2) + draw_grob(ppchleg, x=0.05, y=.2)

save_plot(tar, p, base_height = baseh*1.125, base_width = 3.75, ncol = 3)