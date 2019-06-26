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

# crossings <- function(x, ycmp, tfrle = rle(ycmp > y), super) {
#   sp <- cumsum(head(tfrle$lengths, -1))
#   y <- rep(0, length(ycmp))
#   if (length(sp)) {
#   	expand <- c(1,rep(2, length(sp)))
#   	expand[length(expand)] <- 1
#   	
#     x1   = x[sp]; x2=x[sp+1]
#     yr1  = y[sp]; yr2=y[sp+1]
#     yc1  = ycmp[sp]; yc2=ycmp[sp+1]
#     # yrm  = (yr2-yr1)/(x2-x1); 
#     ycm = (yc2-yc1)/(x2-x1)
#     # yrb  = yr2 - yrm*x2;
#     ycb = yc2 - ycm*x2
#     xint = -ycb/ycm
#     keeps <- inverse.rle(tfrle)
#     return(data.table(
#       x = c(x[keeps], x[!keeps], rep(xint,each=2)),
#       y = c(ycmp[keeps], rep(0, length(x[!keeps])), rep(0,sum(expand))),
#       yref = c(super[keeps], rep(0, length(x[!keeps])), rep(0,sum(expand))),
#       side = c(keeps[keeps], rep(FALSE, length(x[!keeps])), rep(tfrle$values, expand))
#     )[order(x)])
#   } else {
#     return(data.table(x=x[integer()], side=logical()))
#   }
# }

# ribs <- cmb.eff[,{
# 	loregions <- crossings(as.numeric(year), lo, tfrle = rle(lo < 0), pmin(0, hi))
# 	redregions <- loregions[side==TRUE][,.(
# 		year = x, ymin=y, ymax=yref, col="under"
# 	)]
# 	
# 	# where lo regions TRUE - lb is in red region
# 	# ... FALSE - lb is in blue region
# 	
# 	hiregions <- crossings(as.numeric(year), hi, tfrle = rle(hi > 0), pmax(0,lo))
#   redregions <- ribbon_regions(as.numeric(year), lo, tfrle = rle(lo < 0), super=hi)[,.(
#       year = x, ymin=y, ymax=pmin(0,super,na.rm=T), col="under"
#   )]
#   blueregions <- ribbon_regions(as.numeric(year), hi, tfrle = rle(hi > 0), super=lo)[,.(
#     year = x, ymin=pmax(0,super,na.rm=T), ymax=y, col="over"
#   )]
#   rbind(redregions, blueregions)
#   },
#   keyby = .(vaccine, catchup, vc_coverage, measure)
# ]
# ribs <- cmb.eff[order(year),{
#   redregions <- ribbon_regions(as.numeric(year), 0, lo, tfrle = rle(lo < 0), super=hi)[,.(
#       year = x, ymin=y, ymax=pmin(0,super,na.rm=T), col="under"
#   )]
#   blueregions <- ribbon_regions(as.numeric(year), 0, hi, tfrle = rle(hi > 0), super=lo)[,.(
#     year = x, ymin=pmax(0,super,na.rm=T), ymax=y, col="over"
#   )]
#   rbind(redregions, blueregions)
#   },
#   keyby = .(vaccine, catchup, vc_coverage, measure)
# ]

p <- ggplot(cmb.eff) + aes(
  shape=vaccine, color=scenario, size=factor(vc_coverage),
  x=year+1, y=value, group=interaction(vaccine, vc_coverage, catchup, measure)
) + theme_minimal() +
  facet_grid(vaccine + catchup ~ vc_coverage, labeller = facet_labels) +
#  geom_ribbon(aes(color=NULL, ymin=lo, ymax=hi, y=NULL), fill="black", alpha=0.5) +
	geom_rect(
		aes(xmin=year+0.5, xmax=year+1.5, ymin=pmax(0,lo), ymax=hi, color=NULL),
		data=cmb.eff[hi > 0],
		fill = "blue", alpha = 0.5
	) +
	geom_rect(
		aes(xmin=year+0.5, xmax=year+1.5, ymin=lo, ymax=pmin(0,hi), color=NULL),
		data=cmb.eff[lo < 0],
		fill = "red", alpha = 0.5
	) +
  # geom_ribbon(
  #   mapping = aes(fill=col, color=NULL, ymin=ymin, ymax=ymax, y=NULL),
  #   data = ribs[col=="over"],
  #   alpha=0.5, show.legend = F
  # ) +
  # geom_ribbon(
  #   mapping = aes(fill=col, color=NULL, ymin=ymin, ymax=ymax, y=NULL),
  #   data = ribs[col=="under"],
  #   alpha=0.5, show.legend = F
  # ) +
  geom_line(alpha=1, size=vc_sizes["0"]) +
  geom_pchline(dt=cmb.eff[catchup == "routine"], fill="white", alpha=1, show.legend = F) +
  geom_pchline(dt=cmb.eff[catchup != "routine"], fill="black", alpha=1, show.legend = F) +
  scale_year() + scale_y_continuous(name="Interaction") +
  scale_shape_vaccine(guide = "none") +
  scale_color_scenario(guide = "none", value="black") +
  scale_fill_interaction() +
  scale_size_vectorcontrol(guide="none") +
  coord_cartesian(xlim=c(0,40), clip="off") +
  TIRSfacettitle +
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