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
plot.dt <- rbind(
  rbind(vec.eff, vac.eff)[, intervention := factor("single", levels=c("single","combined"), ordered = T) ],
  cmb.eff[, intervention := factor("combined", levels=c("single","combined"), ordered = T) ]
)

# if (grepl("alt", tar)) plot.dt <- rbind(plot.dt, cmb.eff)

limits.dt <- plot.dt[,
  .(value=c(0, 1), year=-1),
  by=scenario
]

scn2_lvls <- with(rbind(expand.grid(scn_lvls[-(1:2)], vac_lvls[2:1]),expand.grid(scn_lvls[2], vac_lvls[3])), paste(Var1, Var2, sep="."))
scn2_labels <- c(with(expand.grid(cu_labels[1:2],vac_labels[2:1]), paste(Var2, Var1, sep=", ")), "75% Vector Control")
scn2_cols <- scn_cols[gsub("^(.+)\\..+$","\\1",scn2_lvls)]
scn2_pchs <- vac_nofill_pchs[gsub("^.+\\.(.+)$","\\1",scn2_lvls)]
names(scn2_labels) <- names(scn2_cols) <- names(scn2_pchs) <- scn2_lvls

scale_color_scenario2 <- scale_generator(
	"color", "Single Intervention", scn2_labels, scn2_cols
)

scale_shape_scenario2 <- scale_generator(
  "shape", "Single Intervention", scn2_labels, scn2_pchs
)

naive.eff.cmdvi <- naive.eff[vaccine == "cmdvi"]
naive.eff.edv <- naive.eff[vaccine == "edv"]
naive.eff.cmdvi.thin <- naive.eff.cmdvi[(((year+1) %% 5 == 0) | year == 0)]
naive.eff.edv.thin <- naive.eff.edv[(((year+1) %% 5 == 0) | year == 0)]

annos <- list(
	annotate("line", x=naive.eff.cmdvi$year+1, y=naive.eff.cmdvi$value, size=vc_sizes["75"], linejoin = "mitre", lineend = "butt", color = "#AAAAFF"),
	annotate("line", x=naive.eff.edv$year+1, y=naive.eff.edv$value, size=vc_sizes["75"], linejoin = "mitre", lineend = "butt", color = "#AAAAFF"),
	annotate("point", x=naive.eff.cmdvi.thin$year+1, y=naive.eff.cmdvi.thin$value, shape=vac_pchs["cmdvi"], size=2, color="#55CC55", fill="#55CC55"),
	annotate("point", x=naive.eff.edv.thin$year+1, y=naive.eff.edv.thin$value, shape=vac_pchs["edv"], size=2, color="#55CC55", fill="#55CC55")
)

# illustrate combined effectiveness
# with coverage, catchup example

p1shared <- ggplot(
  plot.dt[intervention == "single" ]
) + aes(
  x=year+1, y=value, color=interaction(scenario, vaccine),
  shape=interaction(scenario, vaccine), size=factor(vc_coverage),
  group = interaction(scenario, catchup, vaccine, vc_coverage, estimate)
) +
  geom_line(linejoin = "mitre", lineend = "butt") +
  geom_point(data=plot.dt[intervention == "single"][((year+1) %% 5 == 0) | (year == 0)], size=2) +
  scale_size_vectorcontrol(guide = "none") +
  scale_fill_catchup(guide="none", na.value=NA) + theme(
    legend.margin = margin(), legend.spacing = unit(25, "pt"),
    legend.spacing.x = unit(-2,"pt"),
    legend.text = element_text(size=rel(0.4), margin = margin(l=unit(6,"pt"))),
    legend.title = element_text(size=rel(0.5)), legend.title.align = 0.5,
    legend.key.height = unit(1,"pt"),
    legend.box.spacing = unit(2.5, "pt")
  )

p1lines <- p1shared +
  scale_color_scenario2(guide=guide_legend(
    override.aes = list(size=vc_sizes[c("0","0","75")], shape=NA)
  )) +
  scale_shape_scenario2(guide="none")

p1shapes <- p1shared +
  scale_color_scenario2(guide="none") +
  scale_shape_scenario2(guide=guide_legend(
    override.aes = list(linetype=0, size=2, color=scn_cols["vac"])
  ))

p1lleg <- get_legend(p1lines)
p1sleg <- get_legend(p1shapes)

# ggdraw(p1+theme(legend.position = "none")) + draw_grob(p1lleg, x=0.5, y=0.4) + draw_grob(p1sleg, x=0.5, y=0.4)
# ggdraw(p1+theme(legend.position = "none")) + draw_grob(p1sleg, x=0.5, y=0.4)

ribbon.dt <- cmb.eff[
  naive.eff[,.(assume.eff=value, intervention=factor("combined", levels=c("single","combined"), ordered = T)), keyby=key(naive.eff)], on=key(naive.eff),
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

plot2.dt <- ribbon.dt[,.(x=year+1, y=assume.eff, ycmp=value), keyby=.(vc_coverage, vaccine, catchup, scenario, intervention)]

illus_labels <- rbind(
  copy(naive.eff[year==27])[, intervention := factor("single",   levels=c("single","combined"), ordered = T)],
  copy(naive.eff[year==27])[, intervention := factor("combined", levels=c("single","combined"), ordered = T)]
)
illus_labels[vaccine == "edv", value := value + 0.1]
illus_labels[vaccine == "cmdvi" & intervention == "single", value := value - 0.072]
illus_labels[vaccine == "cmdvi" & intervention == "combined", value := value - 0.125]
illus_labels[intervention == "single", lab := paste("Naive 75%", vac_labels[vaccine],sep=" + ") ]
illus_labels[intervention == "combined", lab := ifelse(vaccine == "edv", "Amplification", "Interference") ]

resp <- ggplot(
  plot.dt
) + theme_minimal() + aes(
  x=year+1, y=value, color=interaction(scenario, vaccine),
  shape=vaccine, size=factor(vc_coverage),
  group = interaction(scenario, catchup, vaccine, vc_coverage)
) + facet_grid(intervention ~ ., labeller = labeller(intervention=c(single="Single Interventions",combined="Combined Interventions"))) +
  geom_limits(limits.dt) +
  geom_altribbon(plot2.dt, withlines = F) +
  annos +
  geom_line(linejoin = "mitre", lineend = "butt") +
  geom_point(data=plot.dt[((year+1) %% 5 == 0) | (year == 0)], size=2) +
  geom_text(mapping=aes(label=lab, fill=NULL, color=NULL, size=NULL), data=illus_labels, size=2.5) +
  scale_size_vectorcontrol(guide = "none") +
  scale_color_scenario2(guide = "none") +
  scale_shape_vaccine(values=vac_nofill_pchs, guide = "none") +
  scale_fill_interaction(guide="none") +
  scale_year() +
  scale_effectiveness() +
  coord_cartesian(clip = "off") +
  theme(
    panel.spacing.y = unit(15, "pt"), # panel.spacing.x = unit(15, "pt"),
    plot.margin = margin(t = unit(6,"pt"), r = unit(6,"pt")),
    strip.text.y = element_text(angle=90),
    legend.position = "none"
  )

p <- ggdraw(resp) + draw_grob(p1lleg, x=0.55, y=0.45) + draw_grob(p1sleg, x=0.55, y=0.45)

## TODO shrink vertical gap, slightly shrink font, do labels for theoretical combos?,
#  overflow point off end on right

plotutil(p, h=4.5, w=3, tar)