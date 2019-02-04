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

cmb.eff <- effstats.dt[variable == "syn", .(
  value = warnnonunique(med, variable),
  lo=warnnonunique(lo, variable),
  hi=warnnonunique(hi, variable),
  scenario = trans_scnario(1, 1),
  estimate = "simulated"
), keyby=.(
  vaccine = factor(vaccine, rev(levels(vaccine)), ordered = T), catchup, vc_coverage, year,
  measure = variable
)]

ribbon_regions <- function(x, y, ycmp, tfrle = rle(ycmp > y), super) {
  sp <- cumsum(head(tfrle$lengths, -1))
  y <- rep(y, length(ycmp))
  if (length(sp)) {
    x1   = x[sp]; x2=x[sp+1]
    yr1  = y[sp]; yr2=y[sp+1]
    yc1  = ycmp[sp]; yc2=ycmp[sp+1]
    yrm  = (yr2-yr1)/(x2-x1); ycm = (yc2-yc1)/(x2-x1)
    yrb  = yr2 - yrm*x2; ycb = yc2 - ycm*x2
    xint = (yrb-ycb)/(ycm-yrm); yint = yrm*xint + yrb
    keeps <- inverse.rle(tfrle)
    return(data.table(
      x = c(x[keeps], xint),
      y = c(ycmp[keeps], yint),
      super = c(super[keeps], rep(NA, length(xint)))
    )[order(x)])
  } else {
    if (tfrle$value) {
      return(data.table(x=x, y=ycmp, super=super))
    } else {
      return(data.table(x=x[integer()],y=ycmp[integer()], super=super[integer()]))
    }
  }
}

ribs <- cmb.eff[,{
  redregions <- ribbon_regions(as.numeric(year), 0, lo, tfrle = rle(lo < 0), super=hi)[,.(
      year = x, ymin=y, ymax=pmin(0,super,na.rm=T), col="under"
  )]
  blueregions <- ribbon_regions(as.numeric(year), 0, hi, tfrle = rle(hi > 0), super=lo)[,.(
    year = x, ymin=pmax(0,super,na.rm=T), ymax=y, col="over"
  )]
  rbind(redregions, blueregions)
  },
  keyby = .(vaccine, catchup, vc_coverage, measure)
]

p <- ggplot(cmb.eff) + aes(
  shape=vaccine, color=scenario, size=factor(vc_coverage),
  x=year+1, y=value, group=interaction(vaccine, vc_coverage, catchup, measure)
) + theme_minimal() +
  facet_grid(vaccine + catchup ~ vc_coverage, labeller = facet_labels) +
#  geom_ribbon(aes(color=NULL, ymin=lo, ymax=hi, y=NULL), fill="black", alpha=0.5) +
  geom_ribbon(
    mapping = aes(fill=col, color=NULL, ymin=ymin, ymax=ymax, y=NULL),
    data = ribs[col=="over"],
    alpha=0.5, show.legend = F
  ) +
  geom_ribbon(
    mapping = aes(fill=col, color=NULL, ymin=ymin, ymax=ymax, y=NULL),
    data = ribs[col=="under"],
    alpha=0.5, show.legend = F
  ) +
  geom_line(alpha=1, size=vc_sizes["0"]) +
  geom_point(data=cmb.eff[pchstride(year) & catchup == "routine"], fill="white", alpha=1, size=pchsize) +
  geom_point(data=cmb.eff[pchstride(year) & catchup != "routine"], fill="black", alpha=1, size=pchsize) +
  scale_year() + scale_y_continuous(name="Interaction") +
  # scale_fill_interaction(
  #   guide = gds(1, keyheight=unit(12,"pt"), label.position = "right", direction="vertical", override.aes=list(alpha=c(0.4,0.4)))
  # ) +
  scale_shape_vaccine(guide = "none") +
  scale_color_scenario(guide = "none", value="black") +
  scale_fill_interaction() +
  scale_size_vectorcontrol(guide="none") +
  coord_cartesian(xlim=c(0,40), clip="off") +
  theme(
    legend.margin = margin(), legend.spacing = unit(25, "pt"),
    legend.text = element_text(size=rel(0.5)),
    legend.title = element_text(size=rel(0.6)), legend.title.align = 0.5,
    panel.spacing.y = unit(15, "pt"), panel.spacing.x = unit(15, "pt"),
    strip.text.y = element_text(angle=90),
    legend.key.height = unit(1,"pt"),
    legend.box.spacing = unit(2.5, "pt"),
    legend.position = c(100/120, 0.6) # think this looks best, but can comment out to return to margin
  ) +
  scale_alpha_manual(values=c(delta=int_alpha), guide = "none")

plotutil(p, h=7, w=6, tar)