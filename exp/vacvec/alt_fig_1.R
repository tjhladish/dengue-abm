require(data.table)
require(ggplot2)

args <- c("baseline.rds","intervention.rds","effstats.rds","alt_fig_2.png")
args <- commandArgs(trailingOnly = TRUE)

baseline.dt <- readRDS(args[1])
intervention.dt <- readRDS(args[2])
effstats.rds <- readRDS(args[3])

int.keys <- grep("particle", key(intervention.dt), invert = T, value = T)

intervention.dt[vc == 0, vc_coverage := 0]
vac_mechs <- c("cmdvi","traditional","none")

intervention.dt[, vac_mech := factor(vac_mechs[vac_mech + 1], levels = vac_mechs, ordered = T)]
intervention.dt[vac == 0, vac_mech := factor(vac_mechs[3], levels = vac_mechs, ordered = T)]

catchups <- c("none","catchup")
intervention.dt[, catchup := factor(catchups[catchup+1], levels = catchups, ordered = T) ]

base.inc <- baseline.dt[,{
  qs <- quantile(s, probs = c(0.025, .25, .5, .75, .975))
  names(qs) <- c("lo.lo","lo","med","hi","hi.hi")
  as.list(qs)
}, keyby=.(
  vc = rep(0, length(year)), vac = rep(0, length(year)),
  vc_coverage = rep(0, length(year)),
  vac_mech =factor(rep("none", length(year)), levels = vac_mechs, ordered = T),
  catchup = factor(rep("none", length(year)), levels = catchups, ordered = T),
  year
)]

inte.inc <- intervention.dt[xor(vc,vac), {
  qs <- quantile(s, probs = c(0.025, .25, .5, .75, .975))
  names(qs) <- c("lo.lo","lo","med","hi","hi.hi")
  as.list(qs)  
}, keyby=int.keys ]

scn.lvls <- c("vaccine","vc")

inc.plot.dt <- rbind(
  copy(base.inc)[, scenario := factor(scn.lvls[1], levels=scn.lvls, ordered = T ) ],
  copy(base.inc)[, scenario := factor(scn.lvls[2], levels=scn.lvls, ordered = T ) ],
  inte.inc[, scenario := factor(scn.lvls[vc+1], levels=scn.lvls, ordered = T) ]
)

yucpop <- 18.17734 # 100ks

plot.dt <- inc.plot.dt[, value := med/yucpop ][, measure := "incidence" ][, .(value), by=.(vc_coverage, vac_mech, catchup, year, scenario, measure)]

vac.eff <- effstats.rds[variable == "vac.eff", .(
  value=unique(med), measure="effectiveness", vc_coverage=0, scenario="vaccine"),
  keyby=.(vac_mech = vac_mechs[vac_mech+1], catchup=catchups[catchup+1], year)
]
vec.eff <- effstats.rds[variable == "vec.eff", .(
  value=unique(med), measure="effectiveness", scenario="vc", vac_mech=vac_mechs[3], catchup=catchups[1]),
  keyby=.(vc_coverage, year)
]

twoplot <- rbind(plot.dt, vec.eff, vac.eff)

extend <- copy(twoplot[year == 39])[, year := 40 ]

gds <- function(ord) guide_legend(
  title.position = "top", direction = "horizontal", order = ord,
  label.position = "top"
)

facet_labels <- labeller(
  measure = c(effectiveness="Annual Effectiveness", incidence="Annual Incidence per 100k"),
  scenario = c(vaccine="Vaccine Only",vc="Vector Control Only")
)

sizebase <- 0.5
sizestep <- 0.2
sizes <- seq(from=sizebase, by=0.2, length.out = 4)
names(sizes) <- c(0,25,50,75)

limits.dt <- twoplot[,.(value=c(0, ceiling(max(value))), year=-1), by=measure]

limits.dt[measure == "incidence", value := { ceiling(value/250)*250 }]

p <- ggplot(
#  rbind(twoplot, extend) # for step / segment versions
  twoplot
) + aes(
  # x=year,
  x=year + 1,
  y=value, color=vac_mech, linetype=catchup,
    size=factor(vc_coverage), group=interaction(vc_coverage, catchup, vac_mech)) +
  facet_grid(measure ~ scenario, scales = "free_y", switch = "y", labeller = facet_labels) +
#  geom_segment(mapping = aes(yend=value, xend=year+1)) +
# geom_step() +
  geom_blank(mapping=aes(color=NULL, linetype=NULL, size=NULL, group=NULL), data=limits.dt) +
  geom_line() +
  theme_minimal() +
  scale_color_manual("Vaccine Mechanism", values=c(none="black",cmdvi="blue",traditional="green"), guide=gds(1)) +
  scale_linetype_manual("Catchup", values=c(none="solid", catchup="dashed"), guide=gds(2)) +
  scale_size_manual("Vector Control Coverage %",values=sizes, guide=gds(3)) +
  scale_x_continuous("Year", expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(xlim=c(0,40)) +
  theme(
    axis.title.y = element_blank(),
    strip.placement = "outside",
    legend.box = "horizontal",
    legend.position = c(0.5, 0.5),
    legend.justification = c(0.5, 0.5),
    legend.margin = margin(), legend.spacing = unit(25, "pt"),
    legend.text = element_text(size=rel(0.5)),
    legend.title = element_text(size=rel(0.6)), legend.title.align = 0.5,
    panel.spacing.y = unit(30, "pt"), panel.spacing.x = unit(15, "pt"),
    legend.key.height = unit(1,"pt")
  )

ggsave(
  tail(args,1), p, device = "png",
  width = 6, height = 4, dpi = "retina", units = "in"
)