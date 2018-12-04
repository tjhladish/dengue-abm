require(data.table)
require(ggplot2)

args <- c("effstats.rds")
args <- commandArgs(trailingOnly = TRUE)

stat.eff.dt <- readRDS(args[1])

vac.only <- stat.eff.dt[variable == "vac.eff", .(vac.eff=unique(med)), keyby=.(vac_mech, catchup, year)]
naive <- stat.eff.dt[variable == "ind.eff", .(assume.eff = med), keyby=.(vc_coverage, vac_mech, catchup, year)]

ribbon.dt <- naive[vac.only, on=.(vac_mech, catchup, year), nomatch=0, allow.cartesian=TRUE]

combo.dt <- stat.eff.dt[variable == "combo.eff",.(eff = med), keyby=.(vc_coverage, vac_mech, catchup, year)]

other.ribbon <- combo.dt[naive, on=.(vc_coverage, vac_mech, catchup, year), nomatch=0][vac.only, on=.(vac_mech, catchup, year), nomatch=0, allow.cartesian=T]

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

geom_altribbon <- function(dt, withlines = TRUE, by=NULL) {
  res <- dt[, with(ribbon_intercepts(x, y, ycmp), {
    xlims <- c(x[1], xint, x[.N])
    inner.dt <- rbind(.SD, data.table(x=xint, y=yint, ycmp=yint))[order(x)]
    byvar <- .BY
    inner.dt[, names(byvar) := byvar ]
    # browser()
    res <- lapply(1:(length(xlims)-1), function(i) {
      slice <- inner.dt[between(x, xlims[i], xlims[i+1])]
      geom_polygon(
        aes(x=x, y=y, fill=slice[,any(ycmp>y)], linetype=NA),
        slice[,.(x=c(x,rev(x)), y=c(y,rev(ycmp))), by=names(byvar)],
        alpha = 0.2
      )
    })
    .(polys=res)
  }), by=by ]$polys
  if (withlines) {
    res <- c(res,
      geom_line(aes(x=x, y=y, alpha=NULL, size="assumed"), alpha=0.5, data=dt),
      geom_line(aes(x=x, y=ycmp, alpha=NULL, size="simulated"), data=dt)
      #dt[,.(polys=list(geom_line(aes(x=x, y=ycmp, linetype="comparison"), data=.SD))), by=by]$polys
    )
  }
  res
}

vac_mechs <- c("5th Sero", "Traditional")
vac_cols <- c("green", "blue")
names(vac_cols) <- vac_mechs
vacfact <- function(i) factor(vac_mechs[i+1], levels = rev(vac_mechs), ordered = T)

vac.only[, vc_coverage := 0 ]

p <- ggplot(other.ribbon) + theme_minimal() + aes(group=factor(vc_coverage), linetype=factor(vc_coverage)) +
  geom_altribbon(other.ribbon[,.(x=year, y=assume.eff, ycmp=eff), by=.(vc_coverage, vac_mech, catchup)], by=c("vc_coverage","vac_mech","catchup")) +
  geom_line(aes(x=year, y=vac.eff, size="simulated"), vac.only) +
  facet_grid(vacfact(vac_mech) ~ c("No Catchup","Catchup")[catchup+1]) +
  scale_size_manual("Combination",
    values=c(simulated=0.75, assumed=0.2),
    guide=guide_legend(order=1, direction = "horizontal", title.position = "top")
  ) +
  scale_fill_manual("Interaction",
    labels=c(`FALSE`="inhibit",`TRUE`="enhance"),
    values=c(`FALSE`="red",`TRUE`="blue"),
    guide=guide_legend(direction = "horizontal", title.position = "top")
  ) +
  scale_x_continuous("Year", expand = c(0,0)) + scale_y_continuous("Annual Effectiveness", expand=c(0,0)) +
  coord_cartesian(ylim=c(0,1)) +
  scale_linetype_discrete("Vector Control Coverage",
    guide=guide_legend(
      override.aes = list(fill=NA), order=2,
      direction = "horizontal", title.position = "top"
    )
  ) + theme(
    legend.box = "horizontal",
    legend.position = c(0.50,0.475), legend.justification = c(0.5, 0.5),
    panel.spacing.y = unit(20,"pt")
  )

ggsave(
  tail(args,1), p, device = "png",
  width = 7.5, height = 5, dpi = "retina", units = "in"
)