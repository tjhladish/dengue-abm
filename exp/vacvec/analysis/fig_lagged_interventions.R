suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
})

args <- c("figref.rda", "rds/lag_effstats.rds", "rds/effstats.rds", "fig/fig_5.png")
args <- commandArgs(trailingOnly = TRUE)

load(args[1])
stat.eff.dt <- readRDS(args[2])
ref.stat.dt <- readRDS(args[3])
tar <- args[4]

ekeys <- key(stat.eff.dt)

vac.only <- stat.eff.dt[variable == "vac.eff", .(vac.eff=unique(med)), keyby=.(vaccine, catchup, year, vac_first)]
cvac.only <- stat.eff.dt[variable == "c.vac.eff", .(vac.eff=unique(med)), keyby=.(vaccine, catchup, year, vac_first)]
naive <- stat.eff.dt[variable == "ind.eff", .(assume.eff = med), keyby=.(vac_first, vc_coverage, vaccine, catchup, year)]
cnaive <- stat.eff.dt[variable == "c.ind.eff", .(assume.eff = med), keyby=.(vac_first, vc_coverage, vaccine, catchup, year)]

combo.dt <- stat.eff.dt[variable == "combo.eff",.(eff = med), keyby=.(vac_first, vc_coverage, vaccine, catchup, year)]
ccombo.dt <- stat.eff.dt[variable == "c.combo.eff",.(eff = med), keyby=.(vac_first, vc_coverage, vaccine, catchup, year)]

other.ribbon <- combo.dt[naive,
  on=.(vc_coverage, vaccine, catchup, year, vac_first), nomatch=0
][vac.only,
  on=.(vaccine, catchup, year, vac_first), nomatch=0, allow.cartesian=T
][, measure := factor("annual", levels = c("annual","cumulative"), ordered = T) ]

cother.ribbon <- ccombo.dt[cnaive,
  on=.(vc_coverage, vaccine, catchup, year, vac_first), nomatch=0
][cvac.only,
  on=.(vaccine, catchup, year, vac_first), nomatch=0, allow.cartesian=T
][, measure := factor("cumulative", levels = c("annual","cumulative"), ordered = T) ]

allribbon <- rbind(other.ribbon, cother.ribbon)

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
  # res[[1]]$show.legend <- T
  if (withlines) {
    res <- c(res,
             geom_line(aes(x=x, y=y,    alpha="reference"), data=dt),
             geom_line(aes(x=x, y=ycmp, alpha="compareto"), data=dt)
    )
  }
  res
}

ref.combo <- rbind(
  copy(ref.stat.dt)[, vac_first := 1 ],
  copy(ref.stat.dt)[, vac_first := 0 ]
)[variable %in% c("combo.eff","c.combo.eff") & vaccine == "edv" & catchup == "vc+vac" & vc_coverage == 75]

ref.combo[, measure := ifelse(variable == "combo.eff","annual","cumulative") ]

facet_labs <- labeller(
  vac_first = c(`0`="Vector Control First", `1`="Vaccine First"),
  measure = c(annual="Annual Effectiveness", cumulative="Cumulative Effectiveness")
)

p <- ggplot() + theme_minimal() + #aes(group=vaccine, linetype=factor(vaccine)) +
  geom_altribbon(
    allribbon[,.(x=year+1, y=assume.eff, ycmp=eff), keyby=.(vc_coverage, vaccine, catchup, vac_first, measure)],
    withlines = F
  ) +
  geom_line(mapping=aes(x=year+1, y=med, color="reference"), data=ref.combo) +
  geom_line(mapping=aes(x=year+1, y=eff, color="observed"), data=allribbon) +
  facet_grid_freey(vac_first ~ measure, labeller = facet_labs) +
  scale_fill_interaction(guide="none") +
  scale_year() +
  scale_effectiveness() +
  coord_cartesian(ylim=c(0.5,1), xlim=c(0,40)) +
  theme(
    legend.box = "horizontal",
    legend.position = c(0.5,0.5), legend.justification = c(0.5, 0.5),
    legend.text = element_text(size=rel(0.6)),
    legend.title = element_text(size=rel(0.7), vjust = 0),
    legend.title.align = 0.5,
    panel.spacing.y = unit(30,"pt"),
    panel.spacing.x = unit(15,"pt")
  ) +
  scale_colour_manual(name=NULL,
    values=c(reference="grey",observed="black"),
    labels=c(reference="Simultaneous Reference", observed="Lagged Result")
  ) +
  scale_alpha_manual(values=c(delta=0.2), guide="none")

plotutil(p, h=6, w=6, tar)