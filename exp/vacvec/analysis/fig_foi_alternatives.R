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

ribbon_intercepts <- function(x, y, ycmp) {
	sp <- cumsum(head(rle(ycmp > y)$lengths, -1))
	if (length(sp)) {
		x1   = x[sp]; x2=x[sp+1]
		yr1  = y[sp]; yr2=y[sp+1]
		yc1  = ycmp[sp]; yc2=ycmp[sp+1]
		yrm  = (yr2-yr1)/(x2-x1); ycm = (yc2-yc1)/(x2-x1)
		yrb  = yr2 - yrm*x2; ycb = yc2 - ycm*x2
		xint = (yrb-ycb)/(ycm-yrm); yint = yrm*xint + yrb
		return(list(xint=xint, yint=yint))
	} else return(list(xint=double(),yint=double()))
}

geom_altribbon <- function(dt, ky=key(dt)) {
	res <- dt[, with(ribbon_intercepts(x, y, ycmp), {
		xlims <- c(x[1], xint, x[.N])
		inner.dt <- cbind(rbind(
			.SD,
			data.table(x=xint, y=yint, ycmp=yint)
		)[order(x)], as.data.table(.BY))
		res <- lapply(1:(length(xlims)-1), function(i) {
			slice <- inner.dt[between(x, xlims[i], xlims[i+1])]
			geom_polygon(
				mapping=aes(
					x=x, y=y, linetype=NULL, shape=NULL, color=NULL, size=NULL, alpha="delta",
					fill=col
				),
				data=cbind(slice[,.(
					x=c(x,rev(x)), y=c(y,rev(ycmp)), col=trans_int(slice[,any(ycmp>y)])
				)], as.data.table(.BY)),
				show.legend = F
			)
		})
		.(polys=res)
	}), keyby=ky ]$polys
	res
}

ribbons <- dcast.data.table(stat.eff.dt[variable %in% c("combo.eff","ind.eff")], foi + vaccine + catchup + year + vc_coverage ~ variable, value.var = "med")[,
	.(x=year+1, y=ind.eff, ycmp=combo.eff), keyby=.(foi, vc_coverage, vaccine, catchup)
]

pbase <- ggplot(
  real.dt
) + theme_minimal() +
  aes(shape=vaccine, color=scenario, x=year+1, y=med, size=factor(vc_coverage)) +

	geom_altribbon(ribbons) +
	
  geom_line(mapping=aes(color="vc+naive"), data=stat.eff.dt[variable == "ind.eff"]) +
  geom_pchline(dt=stat.eff.dt[variable == "ind.eff"], color=light_cols["vac"], fill="white") +
  geom_line(data=real.dt[scenario == "vc"]) +

  geom_line(data=real.dt[scenario != "vc"]) +
  geom_pchline(dt=real.dt, fill="white") +
  scale_size_vectorcontrol(guide="none") +
  scale_shape_vaccine(guide="none") +
  scale_year() +
  scale_effectiveness() +
	scale_fill_interaction(guide="none") +
	scale_alpha_manual(values=c(delta=int_alpha), guide = "none") +
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
  	legend.key.height = unit(12,"pt"),
  	legend.box.spacing = unit(2.5, "pt")
  )

labs <- scn_labels
labs["vc"] <- paste0(vc_labels["75"]," ",scn_labels["vc"])
labs["vac"] <- paste0("Routine ", vac_labels["cmdvi"]," Only")
labs["vc+vac"] <- paste0("TIRS"," & ",vac_labels["cmdvi"])
labs["vc+naive"] <- labs["vac+naive"] <- "Naive Estimate"

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

p <- ggdraw(pbase + scale_color_scenario(guide="none")) + draw_grob(pltyleg, x=0.42, y=.18) + draw_grob(ppchleg, x=0.42, y=.18)

save_plot(tar, p, base_height = baseh*1.125, base_width = 3.75, ncol = 3)