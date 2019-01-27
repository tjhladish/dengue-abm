suppressPackageStartupMessages({
  require(data.table)
  require(cowplot)
})

warnnonunique <- function(var, variable, collapse = median) {
  if (length(unique(var)) != 1) warning(sprintf("non unique %s", unique(variable)))
  collapse(var)
}

# debugging args for interactive use
args <- c("figref.rda", "rds/effstats.rds", "fig/fig_3.png")

# expected args:
#  1-3 required: reference_results, interventions_results, effectiveness_stats
#  optional: slice of plot facets
#  required: target plot file
args <- commandArgs(trailingOnly = TRUE)

# load the reference digests
load(args[1])
effstats.dt    <- readRDS(args[2])
tar <- tail(args, 1)

vac.eff <- effstats.dt[variable == "vac.eff" & vc_coverage == 75 & catchup=="vc+vac", .(
    value = warnnonunique(med, variable),
    vc_coverage = 0,
    scenario = trans_scnario(0, 1),
    catchup = "vac-only",
    estimate = "simulated"
  ), keyby=.(
    vaccine, year
    #, measure = trans_meas(gsub("vac\\.","",variable))
  )
]

vec.eff <- cbind(effstats.dt[variable == "vec.eff" & vc_coverage == 75 & catchup == "vc+vac", .(
    value = warnnonunique(med, variable),
    scenario = trans_scnario(1, 0),
    estimate = "simulated"
  ), keyby=.(
    vc_coverage, year
    #, measure = trans_meas(gsub("vec\\.","",variable))
  )
], reference.scenario[,.(vaccine, catchup)])

cmb.eff <- effstats.dt[variable == "combo.eff" & vc_coverage == 75 & catchup == "vc+vac", .(
	value = warnnonunique(med, variable),
	scenario = trans_scnario(1, 1),
	estimate = "simulated"
), keyby=.(vaccine, catchup, vc_coverage, year)]

naive.eff <- effstats.dt[variable == "ind.eff" & vc_coverage == 75 & catchup == "vc+vac", .(
	value = warnnonunique(med, variable),
	scenario = trans_scnario(1, 1),
	estimate = "naive"
), keyby=.(vaccine, catchup, vc_coverage, year)]

# plot.dt <- rbind(vec.eff, vac.eff, cmb.eff)
plot.dt <- rbind(vec.eff, vac.eff, naive.eff)

# if (grepl("alt", tar)) plot.dt <- rbind(plot.dt, cmb.eff)

limits.dt <- plot.dt[,
  .(value=c(0, 1), year=-1),
  by=scenario
]

# illustrate combined effectiveness
# with coverage, catchup example
p1 <- ggplot(
  plot.dt
) + theme_minimal() + aes(
  x=year + 1, y=value, color=scenario,
  fill=catchup, shape=vaccine, size=factor(vc_coverage), alpha = estimate,
  group = interaction(scenario, catchup, vaccine, vc_coverage, estimate)
#  , alpha = estimate
) +
#  facet_grid_freey(scenario ~ ., labeller = facet_labels) +
#  geom_segment(mapping = aes(yend=value, xend=year+1)) +
# geom_step() +
  geom_limits(limits.dt) +
  geom_line(data=plot.dt[estimate != "naive"], linejoin = "mitre", lineend = "butt") +
	geom_line(data=plot.dt[estimate == "naive"], linejoin = "mitre", lineend = "butt", color = scn_cols["vc"]) +
#	geom_line(data=plot.dt[estimate == "naive"], linejoin = "mitre", lineend = "butt", color = scn_cols["vac"], size = vc_sizes["0"]) +
#	geom_point(size=.75) +
	geom_point(data=plot.dt[(((year+1) %% 5 == 0) | year == 0) & estimate != "naive"], size=2) +
	geom_point(data=plot.dt[(((year+1) %% 5 == 0) | year == 0) & estimate == "naive"], size=2, color=scn_cols["vac"], fill=scn_cols["vac"]) +
  scale_size_vectorcontrol(breaks=vc_lvls[c(1,4)], guide="none") +
  scale_color_scenario(guide = "none") +
  scale_shape_vaccine(guide = "none") +
  scale_fill_catchup(guide="none", na.value=NA) +
	scale_alpha_manual(values = c(naive=0.4,simulated=1), guide="none") +
  scale_year() +
  scale_y_continuous("Effectiveness",expand = c(0,0)) +
  theme(
    panel.spacing.y = unit(15, "pt"), # panel.spacing.x = unit(15, "pt"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    plot.margin = margin(r = unit(6,"pt"))
  )

# show interaction for particular example

ribbon.dt <- cmb.eff[
	naive.eff[,.(assume.eff=value), keyby=key(naive.eff)], on=key(naive.eff),
	nomatch=0
]

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

geom_altribbon <- function(dt, withlines = TRUE, ky=key(dt)) {
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
	res[[1]]$show.legend <- T
	if (withlines) {
		res <- c(res,
						 geom_line(aes(x=x, y=y,    alpha="reference"), data=dt),
						 geom_line(aes(x=x, y=ycmp, alpha="compareto"), data=dt)
		)
	}
	res
}

# dt <- ribbon.dt[,.(x=year+1, y=assume.eff, ycmp=value), keyby=.(vc_coverage, vaccine, catchup)]

plot2.dt <- ribbon.dt[,.(x=year+1, y=assume.eff, ycmp=value), keyby=.(vc_coverage, vaccine, catchup, scenario)]

p2 <- ggplot(plot2.dt) + aes(shape=vaccine, color=scenario, size=factor(vc_coverage), x=x, y=y) + theme_minimal() + #aes(linetype=vaccine, shape=vaccine) +
	geom_altribbon(plot2.dt, withlines = F) +
	geom_line(alpha=0.4, color=scn_cols["vc"]) +
	geom_point(data=plot2.dt[(x %% 5 == 0) | x == 1], size=2, alpha=0.4, fill=scn_cols["vac"], color=scn_cols["vac"]) +
	#geom_point(fill=scn_cols["vc+vac"]) +
	geom_line(aes(y=ycmp)) +
#	geom_point(aes(y=ycmp), size=.75, fill=scn_cols["vc+vac"]) +
	geom_point(aes(y=ycmp), data=plot2.dt[(x %% 5 == 0) | x == 1], size=2, fill=scn_cols["vc+vac"]) +
	scale_year() + scale_y_continuous("Effectiveness", expand=c(0,0)) +
	scale_fill_interaction(guide="none"
		#guide = gds(1, keyheight=unit(12,"pt"), label.position = "top", direction="horizontal", override.aes=list(alpha=c(0.4,0.4)))
	) +
	scale_shape_vaccine(guide = "none") +
	scale_color_scenario(guide = "none") +
	scale_size_vectorcontrol(guide="none") +
	coord_cartesian(ylim=c(0,1), xlim=c(0,40)) +
	theme(
		legend.margin = margin(), legend.spacing = unit(25, "pt"),
		legend.text = element_text(size=rel(0.5)),
		legend.title = element_text(size=rel(0.6)), legend.title.align = 0.5,
		#panel.spacing.y = unit(15, "pt"), # panel.spacing.x = unit(15, "pt"),
		legend.key.height = unit(1,"pt"),
		legend.box.spacing = unit(2.5, "pt"),
		legend.position = c(0.5, .1),
		plot.margin = margin(t = unit(6,"pt"), r = unit(6,"pt"))
	) +
	scale_alpha_manual(values=c(delta=0.4), guide = "none")

## TODO shrink vertical gap, slightly shrink font, do labels for theoretical combos?,
#  overflow point off end on right

plotutil(plot_grid(p1, p2, align = "hv", nrow = 2), h=4.5, w=3, tar)