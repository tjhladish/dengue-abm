suppressPackageStartupMessages({
	require(data.table)
  require(scales)
  require(ggplot2)
	require(cowplot)
}) 

warnnonunique <- function(var, variable, collapse = median) {
	if (length(unique(var)) != 1) warning(sprintf("non unique %s", unique(variable)))
	collapse(var)
}

args <- c("figref.rda", "rds/alt_effstats.rds","fig/fig_4_alt.png")
args <- commandArgs(trailingOnly = TRUE)

load(args[1])
effstats.dt    <- readRDS(args[2])#[eval(mainfilter)]

effstats.dt[, vaccine := factor(as.character(vaccine), levels = c("d50e","d70e","d90e"), ordered = T)]

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

vac.eff <- effstats.dt[variable == "vac.eff", .(
  value = warnnonunique(med, variable),
  vc_coverage = 0,
  scenario = trans_scnario(0, 1)
), keyby=.(
  vaccine, catchup = ifelse(catchup=="vc+vac","vac-only","routine"), year
  #, measure = trans_meas(gsub("vac\\.","",variable))
)
]

vec.eff <- cbind(effstats.dt[variable == "vec.eff", .(
  value = warnnonunique(med, variable),
  scenario = trans_scnario(1, 0)
), keyby=.(
  vc_coverage, year
  #, measure = trans_meas(gsub("vec\\.","",variable))
)
], reference.scenario[,.(vaccine, catchup)])

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

leg.sz <- 0.9

shared <- list(
	facet_grid(vaccine ~ vc_coverage, labeller = labeller(
	  vc_coverage = vc_labels,
	  vaccine = function(ve) paste0(gsub("d(\\d{2})e","\\1%",ve),c("\n","\nDurable Vaccine Efficacy","\n"))
	)),
		geom_line(aes(y=ycmp), alpha=1, size=vc_sizes["0"]),
		geom_pchline(dt=plot2.dt[catchup == "routine"], mapping=aes(y=ycmp), var=expression(x-1), fill="white", alpha=1),
		geom_pchline(dt=plot2.dt[catchup != "routine"], mapping=aes(y=ycmp), var=expression(x-1), fill="black", alpha=1),
		scale_year(), scale_effectiveness(),
		scale_color_scenario(guide = "none", value=c("black")),
		scale_size_vectorcontrol(guide="none"),
		coord_cartesian(ylim=c(0,1), xlim=c(0,40), clip="off"),
		theme_minimal(),
		TIRSfacettitle,
		theme(
			axis.title = element_text(size=rel(1.1)),
			axis.text = element_text(size=rel(1)),
			legend.margin = margin(), legend.spacing = unit(25, "pt"),
			legend.text = element_text(size=rel(leg.sz)),
			legend.title = element_text(size=rel(leg.sz)), legend.title.align = 0.5,
			panel.spacing.y = unit(15, "pt"), panel.spacing.x = unit(15, "pt"),
			strip.text = element_text(size=rel(1.1)),
			strip.text.y = element_text(angle=90),
			legend.key.height = unit(1,"pt"),
			legend.box.spacing = unit(2.5, "pt")
		),
		scale_alpha_manual(values=c(delta=int_alpha), guide = "none")
)

pbase <- ggplot(plot2.dt) + aes(
	shape="d70e",
	color=scenario, #size=factor(vc_coverage),
	x=x, y=y, group=interaction(vaccine, vc_coverage, catchup)
)

# relabs <- cu_labels[c("vc+vac","routine")]
# names(relabs) <- c("d70e","t+cydtdv")
# # get legends
# d70e.lab <- get_legend(pbase + shared + scale_shape_vaccine(vac_labels["d70e"], labels=relabs, guide=guide_legend(
# 	override.aes = list(shape = vac_pchs["d70e"], fill=c(scn_cols["vc+vac"], "white"))
# )))
# cydtdv.lab <- get_legend(pbase + shared + scale_shape_vaccine(vac_labels["t+cydtdv"], labels = relabs, guide=guide_legend(
# 	override.aes = list(shape = vac_pchs["t+cydtdv"], fill=c(scn_cols["vc+vac"], "white"))
# )))

lum <- 60

pmost <- pbase +
  geom_line(
    aes(x=year, y=value, group=vc_coverage),
    data = vec.eff[,.SD, .SDcols=-c("vaccine")],
    color = scales::muted(scn_cols["vc"], lum)
  ) +
  geom_line(
    aes(x=year, y=value, group=interaction(vaccine, catchup)),
    data = vac.eff[,.SD, .SDcols=-c("vc_coverage")],
    color = scales::muted(scn_cols["vac"], lum)
  ) +
  geom_pchline(
    dt = vac.eff[catchup == "routine",.SD, .SDcols=-c("vc_coverage")],
    mapping=aes(x=year, y=value, group=interaction(vaccine, catchup)),
    var=expression(year-1),
    color = scales::muted(scn_cols["vac"], lum),
    fill="white", alpha=1
  ) +
  geom_pchline(
    dt = vac.eff[catchup != "routine",.SD, .SDcols=-c("vc_coverage")],
    mapping = aes(x=year, y=value, group=interaction(vaccine, catchup)),
    var=expression(year-1),
    color = scales::muted(scn_cols["vac"], lum),
    fill = scales::muted(scn_cols["vac"], lum), alpha=1
  ) +
  geom_altribbon(plot2.dt, withlines = F) + shared + scale_shape_vaccine(guide="none") + 
  scale_fill_interaction(guide = "none") + theme(
  legend.position = c(100/120*.285, 0.40) # think this looks best, but can comment out to return to margin
)

save_plot(tar, pmost, ncol = 3, nrow = 3, base_width = 4, base_height = baseh*1.5)
