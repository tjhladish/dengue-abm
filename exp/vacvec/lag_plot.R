require(data.table)
require(ggplot2)

source("projref.R")

args <- c("lag_effstats.rds", "effstats.rds")
args <- commandArgs(trailingOnly = TRUE)

stat.eff.dt <- readRDS(args[1])
ref.stat.dt <- readRDS(args[2])

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
    res <- c(
      res,
      geom_line(aes(x=x, y=y, size="assumed"), alpha=0.5, data=dt),
      geom_line(aes(x=x, y=ycmp, size="simulated"), data=dt)
    )
  }
  res
}

ref.combo <- rbind(
  copy(ref.stat.dt)[, vac_first := 1 ],
  copy(ref.stat.dt)[, vac_first := 0 ]
)[variable %in% c("combo.eff","c.combo.eff") & vaccine == "traditional" & catchup == "catchup" & vc_coverage == 75]

ref.combo[, measure := ifelse(variable == "combo.eff","annual","cumulative") ]

facet_labs <- labeller(
  vac_first = c(`0`="Vector Control First", `1`="Vaccine First"),
  measure = c(annual="Annual", cumulative="Cumulative")
)

gds <- function(order,
                title.position = "top",
                direction = "horizontal",
                ...
) guide_legend(
  title.position = title.position, direction = direction, order = order,
  label.position = "top", ...
)

p <- ggplot(other.ribbon) + theme_minimal() + #aes(group=vaccine, linetype=factor(vaccine)) +
  geom_line(mapping=aes(x=year+1, y=med, size="simulated", color="reference"), data=ref.combo) +
  geom_altribbon(allribbon[,.(x=year+1, y=assume.eff, ycmp=eff), by=.(vc_coverage, vaccine, catchup, vac_first, measure)], by=c("vc_coverage","vaccine","catchup","vac_first","measure")) +
  facet_grid(vac_first ~ measure, labeller = facet_labs) +
  scale_size_manual("Combination",
                    values=c(simulated=1, assumed=0.5),
                    guide=gds(order=1, override.aes=list(fill=NA), keyheight = unit(1,"pt"))
  ) +
  scale_fill_manual("Interaction",
                    labels=c(`FALSE`="interfere",`TRUE`="enhance"),
                    values=c(`FALSE`="red",`TRUE`="blue"),
                    guide=gds(order=3, override.aes=list(alpha=0.2/2), keyheight = unit(5,"pt"))
  ) +
  scale_x_continuous("Year", expand = c(0,0)) +
  scale_y_continuous("Effectiveness", expand=c(0,0)) +
  coord_cartesian(ylim=c(0.5,1), xlim=c(0,40)) +
  theme(
    legend.box = "horizontal",
    legend.position = c(0.5,0.5), legend.justification = c(0.5, 0.5),
    legend.text = element_text(size=rel(0.6)),
    legend.title = element_text(size=rel(0.7), vjust = 0),
    legend.title.align = 0.5,
    panel.spacing.y = unit(30,"pt"),
    panel.spacing.x = unit(15,"pt"),
    strip.text.y = element_text(angle=90)
  ) +
  scale_colour_manual(name=NULL,
    values=c(`reference`="grey"), labels=c(reference="Simultaneous Reference"),
    guide = guide_legend(
      override.aes = list(fill=NA, size=1)
    )
  ) +
  scale_alpha_manual(values=c(delta=0.2), guide="none")

ggsave(
  tail(args,1), p, device = "png",
  width = 7.5, height = 5, dpi = "retina", units = "in"
)