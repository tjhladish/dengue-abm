require(data.table)
require(ggplot2)

args <- c("effstats.rds")
args <- commandArgs(trailingOnly = TRUE)

stat.eff.dt <- readRDS(args[1])

vac.only <- stat.eff.dt[variable == "vac.eff", .(vac.eff=unique(med)), keyby=.(vac_mech, catchup, year)]
naive <- stat.eff.dt[variable == "ind.eff", .(assume.eff = med), keyby=.(vc_coverage, vac_mech, catchup, year)]

ribbon.dt <- naive[vac.only, on=.(vac_mech, catchup, year)]

combo.dt <- stat.eff.dt[variable == "combo.eff",.(eff = med), keyby=.(vc_coverage, vac_mech, catchup, year)]

other.ribbon <- combo.dt[naive, on=.(vc_coverage, vac_mech, catchup, year)][vac.only, on=.(vac_mech, catchup, year)]

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
        aes(x=x, y=y, fill=slice[,any(ycmp>y)], linetype=NA, alpha="delta"),
        slice[,.(x=c(x,rev(x)), y=c(y,rev(ycmp))), by=names(byvar)],
        show.legend = (i==1)
      )
    })
    .(polys=res)
  }), by=by ]$polys
  if (withlines) {
    res <- c(res,
      geom_line(aes(x=x, y=y, size="assumed"), alpha=0.5, data=dt),
      geom_line(aes(x=x, y=ycmp, size="simulated"), data=dt)
      #dt[,.(polys=list(geom_line(aes(x=x, y=ycmp, linetype="comparison"), data=.SD))), by=by]$polys
    )
  }
  res
}

vac_mechs <- c("5th Sero", "Traditional")
vac_cols <- c("green", "blue")
names(vac_cols) <- vac_mechs
vacfact <- function(i) factor(vac_mechs[i+1], levels = rev(vac_mechs), ordered = T)

vac.only.lines <- rbind(
  copy(vac.only)[, vc_coverage := 25 ],
  copy(vac.only)[, vc_coverage := 50 ],
  copy(vac.only)[, vc_coverage := 75 ]
)

facet_labs <- labeller(
  catchup = c(`0`="No Catchup", `1`="Catchup"),
  vc_coverage = c(`25`="\n25%", `50`="Vector Control Coverage\n50%", `75`="\n75%")
)

gds <- function(order,
  title.position = "top",
  direction = "horizontal",
  ...
) guide_legend(
  title.position = title.position, direction = direction, order = order,
  label.position = "top", ...
)

p <- ggplot(other.ribbon) + theme_minimal() + aes(group=vac_mech, linetype=factor(vac_mech)) +
  geom_line(aes(x=year+1, y=vac.eff, size="simulated", color="reference"), vac.only.lines) +
  geom_altribbon(other.ribbon[,.(x=year+1, y=assume.eff, ycmp=eff), by=.(vc_coverage, vac_mech, catchup)], by=c("vc_coverage","vac_mech","catchup")) +
  facet_grid(catchup ~ vc_coverage, labeller = facet_labs) +
  scale_size_manual("Combination",
    values=c(simulated=0.5, assumed=0.2),
    guide=gds(order=1, override.aes=list(fill=NA), keyheight = unit(1,"pt"))
  ) +
  scale_fill_manual("Interaction",
    labels=c(`FALSE`="inhibit",`TRUE`="enhance"),
    values=c(`FALSE`="red",`TRUE`="blue"),
    guide=gds(order=3, override.aes=list(alpha=0.2/6), keyheight = unit(5,"pt"))
  ) +
  scale_x_continuous("Year", expand = c(0,0)) + scale_y_continuous("Annual Effectiveness", expand=c(0,0)) +
  coord_cartesian(ylim=c(0,1), xlim=c(0,40)) +
  scale_linetype_discrete(
    "Vaccine",
    labels=c(`0`="CMDVI",`1`="Traditional"),
    guide=gds(order=2, override.aes = list(fill=NA), keyheight = unit(1,"pt"))
  ) + theme(
    legend.box = "horizontal",
    legend.position = c(0.5,0.5), legend.justification = c(0.5, 0.5),
    panel.spacing.y = unit(30,"pt"),
    panel.spacing.x = unit(15,"pt"),
    strip.text.y = element_text(angle=90),
    legend.text = element_text(size=rel(0.6)),
    legend.title = element_text(size=rel(0.7)),
    legend.title.align = 0.5
    # legend.key = element_rect(),
    # legend.key.height = unit(1, "pt")
#    , strip.background = element_rect(fill="lightgrey", color=NA)
  ) +
  scale_colour_manual(name=NULL,
    values=c(`reference`="grey"), labels=c(reference="Vaccine-only Reference"),
    guide = guide_legend(
      override.aes = list(fill=NA)
    )
  ) +
  scale_alpha_manual(values=c(delta=0.2), guide="none")

ggsave(
  tail(args,1), p, device = "png",
  width = 7.5, height = 5, dpi = "retina", units = "in"
)