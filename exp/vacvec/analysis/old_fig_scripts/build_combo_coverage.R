## building version of combo coverage facet plot

require(data.table)
require(ggplot2)
require(gganimate)

source("projref.R")
source("utils.R")

args <- c("effstats.rds")
args <- commandArgs(trailingOnly = TRUE)

stat.eff.dt <- readRDS(args[1])

vac.only <- stat.eff.dt[variable == "vac.eff", .(eff=median(med)), keyby=.(vaccine, catchup, year)]
vec.only <- stat.eff.dt[variable == "vec.eff", .(eff=median(med)), keyby=.(vc_coverage, year)]
# naive <- stat.eff.dt[variable == "ind.eff", .(assume.eff = med), keyby=.(vc_coverage, vaccine, catchup, year)]
# 
# ribbon.dt <- naive[vac.only, on=.(vaccine, catchup, year), nomatch=0, allow.cartesian=T]
# 
# combo.dt <- stat.eff.dt[variable == "combo.eff",.(eff = med), keyby=.(vc_coverage, vaccine, catchup, year)]
# 
# other.ribbon <- combo.dt[naive, on=.(vc_coverage, vaccine, catchup, year), nomatch=0][vac.only, on=.(vaccine, catchup, year), nomatch=0, allow.cartesian=T]

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

vac.only.lines <- rbind(
  copy(vac.only)[, vc_coverage := 0 ],
  copy(vac.only)[, vc_coverage := 25 ],
  copy(vac.only)[, vc_coverage := 50 ],
  copy(vac.only)[, vc_coverage := 75 ]
)

facet_labs <- labeller(
  catchup = c(none="No Catchup", catchup="Catchup"),
  vc_coverage = c(`0`="No\nVector Control", `25`="\n25%", `50`="Vector Control Coverage\n50%", `75`="\n75%")
)

gds <- function(
  order,
  title.position = "top",
  direction = "horizontal",
  ...
) guide_legend(
  title.position = title.position, direction = direction, order = order,
  label.position = "top", ...
)

consistent <- list(
  theme_minimal(),
  theme(
    legend.box = "horizontal",
    legend.justification = c(0.5, 0.5),
    legend.text = element_text(size=rel(0.6)),
    legend.title = element_text(size=rel(0.7), vjust = 0),
    legend.title.align = 0.5,
    panel.spacing.y = unit(30,"pt"),
    panel.spacing.x = unit(15,"pt"),
    strip.text.y = element_text(angle=90)
    # legend.key = element_rect(),
    # legend.key.height = unit(1, "pt")
    #    , strip.background = element_rect(fill="lightgrey", color=NA)
  ),
  scale_x_continuous("Year", expand = c(0,0)),
  scale_y_continuous("Annual Effectiveness", expand = c(0,0)),
  scale_linetype_manual("Vaccine",
    values=c(cmdvi="44", traditional="solid"),
    guide=gds(order=2, override.aes = list(fill=NA), keyheight = unit(1,"pt"))
  ),
  facet_grid(catchup ~ vc_coverage, labeller = facet_labs),
  coord_cartesian(ylim=c(0,1),xlim=c(0,40))
)

altrenderer <- gifski_renderer(loop = FALSE)

# start with vaccine (both mechs)
vaccine.p <- ggplot(
  vac.only.lines[catchup=="none" & vc_coverage == 0]
) + aes(linetype = vaccine, x=year+1, y=eff) +
  geom_line(mapping = aes(color="foreground")) + geom_point() + 
  geom_rug(sides = "l") +
  geom_text(aes(x = 0, label = sprintf("%.2f",eff)), nudge_y = 0.01, hjust = 0, vjust = 0) +
  consistent + theme(
    legend.position = "bottom"
  ) +
  scale_color_manual(guide = "none", values=c(foreground="black")) +
  transition_reveal(year) +
  ease_aes('linear')

gganimate::anim_save("vaccine_only.gif", vaccine.p, renderer = altrenderer)

catchup.p <- ggplot(
  vac.only.lines[vc_coverage == 0]
) + aes(linetype = vaccine, x=year+1, y=eff) +
  geom_line(mapping = aes(color="foreground")) + geom_point() + 
  geom_rug(sides = "l") +
  geom_text(aes(x = 0, label = sprintf("%.2f",eff)), nudge_y = 0.01, hjust = 0, vjust = 0) +
  consistent + theme(
    legend.position = c(.5,.5)
  ) +
  scale_color_manual(guide = "none", values=c(foreground="black")) +
  transition_reveal(year) +
  ease_aes('linear')

gganimate::anim_save("vaccine_w_catchup.gif", catchup.p, renderer = altrenderer)

vec.p <- ggplot(vec.only[vc_coverage == 50]) +
  aes(linetype = vaccine, x=year+1, y=eff) +
  geom_line(
    data = vac.only.lines[vc_coverage == 50, .(eff), by=.(
      vaccine, t=year+1, catchup, vc_coverage
    )],
    mapping = aes(color="background", x=t)
  ) + 
  geom_line(mapping = aes(color="foreground")) +
  geom_point() + 
  geom_rug(sides = "l") +
  geom_text(aes(x = 0, label = sprintf("%.2f",vac.eff)), nudge_y = 0.01, hjust = 0, vjust = 0) +
  consistent + theme(
    legend.position = c(.5,.5)
  ) +
  scale_color_manual(guide = "none", values=c(foreground="black", background="grey")) +
  transition_reveal(year) +
  ease_aes('linear')

gganimate::anim_save("vector_control.gif", vec.p, renderer = altrenderer)

# then split vertical to show catchup campaigns
# then shadow vector control (50%?)
# then show "expected" combination
# then actual combination
# then difference shading
# then split to 25% / 75%


p <- ggplot(other.ribbon) + theme_minimal() + aes(group=vaccine, linetype=factor(vaccine)) +
  geom_line(aes(x=year+1, y=vac.eff, size="simulated", color="reference"), vac.only.lines) +
  geom_altribbon(other.ribbon[,.(x=year+1, y=assume.eff, ycmp=eff), by=.(vc_coverage, vaccine, catchup)], by=c("vc_coverage","vaccine","catchup")) +
  facet_grid(catchup ~ vc_coverage, labeller = facet_labs) +
  scale_size_manual("Combination",
                    values=c(simulated=1, assumed=0.5),
                    guide=gds(order=1, override.aes=list(fill=NA), keyheight = unit(1,"pt"))
  ) +
  scale_fill_manual("Interaction",
                    labels=c(`FALSE`="interfere",`TRUE`="enhance"),
                    values=c(`FALSE`="red",`TRUE`="blue"),
                    guide=gds(order=3, override.aes=list(alpha=0.2/6), keyheight = unit(5,"pt"))
  ) +
  scale_x_continuous("Year", expand = c(0,0)) + scale_y_continuous("Annual Effectiveness", expand=c(0,0)) +
  coord_cartesian(ylim=c(0,1), xlim=c(0,40)) +
  scale_linetype_discrete(
    "Vaccine",
    labels=c(cmdvi="CMDVI",traditional="Traditional"),
    guide=gds(order=2, override.aes = list(fill=NA), keyheight = unit(1,"pt"))
  ) + theme(
    legend.box = "horizontal",
    legend.position = c(0.5,0.5), legend.justification = c(0.5, 0.5),
    legend.text = element_text(size=rel(0.6)),
    legend.title = element_text(size=rel(0.7), vjust = 0),
    legend.title.align = 0.5,
    panel.spacing.y = unit(30,"pt"),
    panel.spacing.x = unit(15,"pt"),
    strip.text.y = element_text(angle=90)
    # legend.key = element_rect(),
    # legend.key.height = unit(1, "pt")
    #    , strip.background = element_rect(fill="lightgrey", color=NA)
  ) +
  scale_colour_manual(name=NULL,
                      values=c(`reference`="grey"), labels=c(reference="Vaccine-only Reference"),
                      guide = guide_legend(
                        override.aes = list(fill=NA, size=1)
                      )
  ) +
  scale_alpha_manual(values=c(delta=0.2), guide="none")

ggsave(
  tail(args,1), p, device = "png",
  width = 7.5, height = 5, dpi = "retina", units = "in"
)