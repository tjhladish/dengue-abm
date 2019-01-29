suppressPackageStartupMessages({
	require(data.table)
	require(cowplot)
}) 

warnnonunique <- function(var, variable, collapse = median) {
	if (length(unique(var)) != 1) warning(sprintf("non unique %s", unique(variable)))
	collapse(var)
}

args <- c("figref.rda", "rds/effstats.rds","fig/fig_4.png")
args <- commandArgs(trailingOnly = TRUE)

load(args[1])
effstats.dt <- readRDS(args[2])
tar <- tail(args, 1)

cmb.eff <- effstats.dt[variable == "combo.eff", .(
	value = warnnonunique(med, variable),
	scenario = trans_scnario(1, 1),
	estimate = "simulated"
), keyby=.(vaccine, catchup, vc_coverage, year)]

naive.eff <- effstats.dt[variable == "ind.eff", .(
	value = warnnonunique(med, variable),
	scenario = trans_scnario(1, 1),
	estimate = "naive"
), keyby=.(vaccine, catchup, vc_coverage, year)]

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

plot2.dt <- ribbon.dt[,.(x=year+1, y=assume.eff, ycmp=value), keyby=.(vc_coverage, vaccine, catchup, scenario)]

# plot2.dt$vaccine <- factor(plot2.dt$vaccine, levels=rev(levels(plot2.dt$vaccine)), ordered = T)

p <- ggplot(plot2.dt) + aes(
	shape=vaccine, color=scenario, size=factor(vc_coverage),
	x=x, y=y, group=interaction(vaccine, vc_coverage, catchup)
) + theme_minimal() +
	facet_grid(vaccine ~ vc_coverage, labeller = facet_labels) +
	geom_altribbon(plot2.dt, withlines = F) +
	geom_line(aes(y=ycmp), alpha=1, size=vc_sizes["0"]) +
  geom_point(aes(y=ycmp), data=plot2.dt[((x %% 5 == 0) | x == 1) & catchup == "routine"], fill="white", alpha=1, size=2) +
  geom_point(aes(y=ycmp), data=plot2.dt[((x %% 5 == 0) | x == 1) & catchup != "routine"], fill="black", alpha=1, size=2) +
	scale_year() + scale_y_continuous("Effectiveness", expand=c(0,0)) +
	scale_fill_interaction(
		guide = gds(1, keyheight=unit(12,"pt"), label.position = "right", direction="vertical", override.aes=list(alpha=c(0.4,0.4)))
	) +
	scale_pchlty_vaccine(guide = "none") +
	scale_color_scenario(guide = "none", value="black") +
	scale_size_vectorcontrol(guide="none") +
	coord_cartesian(ylim=c(0,1), xlim=c(0,40), clip="off") +
	theme(
		legend.margin = margin(), legend.spacing = unit(25, "pt"),
		legend.text = element_text(size=rel(0.5)),
		legend.title = element_text(size=rel(0.6)), legend.title.align = 0.5,
		panel.spacing.y = unit(15, "pt"), panel.spacing.x = unit(15, "pt"),
		legend.key.height = unit(1,"pt"),
		legend.box.spacing = unit(2.5, "pt")
	) +
	scale_alpha_manual(values=c(delta=int_alpha), guide = "none")

# TODO: switch vaccine precedence

plotutil(p, h=5, w=7.5, tar)