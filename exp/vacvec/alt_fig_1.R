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

gds <- function(ord) guide_legend(title.position = "top", direction = "horizontal", order = ord)

facet_labels <- labeller(
  measure = c(effectiveness="Annual Effectiveness", incidence="Annual Incidence per 100k"),
  scenario = c(vaccine="Vaccine Only",vc="Vector Control Only")
)

p<-ggplot(rbind(twoplot, extend)) + aes(x=year, y=value, color=vac_mech, linetype=catchup,
    size=factor(vc_coverage), group=interaction(vc_coverage, catchup, vac_mech)) +
  facet_grid(measure ~ scenario, scales = "free_y", switch = "y", labeller = facet_labels) +
  geom_step() + theme_minimal() +
  scale_color_manual("Vaccine Mechanism", values=c(none="black",cmdvi="blue",traditional="green"), guide=gds(1)) +
  scale_linetype_manual("Catchup", values=c(none="solid", catchup="dotted"), guide=gds(2)) +
  scale_size_manual("Vector Control Coverage %",values=c(`0`=0.5,`25`=0.75,`50`=1,`75`=1.25), guide=gds(3)) +
  theme(
    axis.title.y = element_blank(),
    strip.placement = "outside",
    legend.box = "horizontal",
    legend.position = c(0.5, 0.5),
    legend.justification = c(0.5, 0.5),
    legend.margin = margin(), legend.spacing = unit(25, "pt"),
    legend.text = element_text(size=rel(0.5)),
    legend.title = element_text(size=rel(0.6)),
    panel.spacing.y = unit(40, "pt")
  )

ggsave(
  tail(args,1), p, device = "png",
  width = 6, height = 4, dpi = "retina", units = "in"
)