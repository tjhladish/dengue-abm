suppressPackageStartupMessages({
  require(data.table)
  require(cowplot)
  source("utils.R")
  source("labelref.R")
})

warnnonunique <- function(var, variable, collapse = median) {
  if (length(unique(var)) != 1) warning(sprintf("non unique %s", unique(variable)))
  collapse(var)
}

args <- c("effstats.rds","fig_2.png")
args <- commandArgs(trailingOnly = TRUE)

stat.eff.dt <- readRDS(args[1])
skeys <- key(stat.eff.dt)

stat.eff.dt[, measure := trans_meas(variable) ]

limits <- cbind(data.table(
  measure = trans_meas(c(rep("combo.eff", 2), rep("syn", 2))),
  year = -1,
  med = c(c(0, 1), c(-.15,.25))
), reference.scenario[,.(catchup, vaccine)])

# vac_only_reference <- stat.eff.dt[variable == "vac.eff",
#   .(measure = trans_meas("combo.eff"),
#     med=warnnonunique(med, variable),
#     vc_coverage = 0),
#   keyby=c(setdiff(skeys, "vc_coverage"))
# ]

legend_theme <- theme(
  legend.margin = margin(), legend.spacing = unit(25, "pt"),
  legend.text = element_text(size=rel(0.5)),
  legend.title = element_text(size=rel(0.6), margin = margin()),
  legend.title.align = 0.5,
  legend.box.margin = margin(),
  panel.spacing.y = unit(15, "pt"),
  legend.key.height = unit(1,"pt"),
  legend.direction = "horizontal", legend.box = "horizontal",
  legend.justification = c(0.5, 1),
  legend.position = c(0.5, -0.1),
  axis.title.x = element_text(size=rel(0.7), margin = margin(t=1, b=30))
)

# perspectives
#  - all combination results
#  - one scenario (75% coverage, catchup, EDV),
#    showing how combination effectiveness compares naive combo of single ints
#  - ibid, but for dengvaxia

allres <- stat.eff.dt[measure == "combo.eff"][,
  perspective := "all"
][,
  type := "obs"
]

synergy.dt <- function(dt, taryear) {
  rel.dt <- dcast.data.table(
    dt[year == taryear], . ~ variable,
    value.var = "med"
  )
  rel.dt <- rbind(rel.dt, rel.dt)
  rel.dt[, arr := c("ref","syn") ]
  rel.dt[ arr == "ref", c("start", "end") := .(max(vec.eff, vac.eff), ind.eff) ]
  rel.dt[ arr == "syn", c("start", "end") := .(ind.eff, combo.eff) ]
  rel.dt[ arr == "syn", arr := ifelse(ind.eff > combo.eff, "syn-","syn+")][,
    .(start, end), by=arr
  ]
  # black arrow from whichever is higher of single curves
  # to grey curve
  # ? plus mark at mid between single curves
  # colored arrow from ind.eff to combo.eff.
  # color blue if positive, red if negative
  # double width of black arrow
}

ty <- 20 # taryear

perspective1 <- stat.eff.dt[
  vaccine == "traditional" &
  vc_coverage == 75 &
  catchup == "catchup" &
  variable %in% c("combo.eff","vac.eff","vec.eff","ind.eff")
][, perspective := "edv" ][, type := "obs" ]
perspective1[variable == "vac.eff", vc_coverage := 0 ]
perspective1[variable == "vec.eff", vaccine := trans_vac(0,0) ]
perspective1[variable == "ind.eff", type := "est" ]
arr1 <- synergy.dt(perspective1, ty-1)[, perspective := "edv" ]

perspective2 <- stat.eff.dt[
  vaccine == "cmdvi" &
  vc_coverage == 75 &
  catchup == "catchup" &
  variable %in% c("combo.eff","vac.eff","vec.eff","ind.eff")
][, perspective := "cmdvi" ][, type := "obs" ]
perspective2[variable == "vac.eff", vc_coverage := 0 ]
perspective2[variable == "vec.eff", vaccine := trans_vac(0,0) ]
perspective2[variable == "ind.eff", type := "est" ]
arr2 <- synergy.dt(perspective2, ty-1)[, perspective := "cmdvi" ]

pers_labels <- labeller(
  perspective = c(
    all="All\nCombination Interventions",
    cmdvi="Dengvaxia-like, Catchup &\n75% Vector Control",
    edv="EDV, Catchup &\n75% Vector Control"
  )
)

p<-ggplot(
  rbind(allres, perspective1, perspective2)
  #stat.eff.dt[measure %in% c("combo.eff")]
) +
  theme_minimal() +
  facet_grid(. ~ perspective, labeller=pers_labels) +
  aes(x=year+1, y=med,
    linetype=vaccine, size=factor(vc_coverage),
    alpha = type, color=catchup,
    group = interaction(variable, vaccine, vc_coverage, catchup, type)
  ) +
  geom_segment(
    data = rbind(arr1, arr2)[arr == "syn+"],
    mapping = aes(x=ty, xend=ty, y=start, yend=end), inherit.aes = FALSE,
    arrow = arrow(15,length = unit(0.05,"npc"), type = "closed"),
    color = "blue", size = 2
  ) +
  geom_segment(
    data = rbind(arr1, arr2)[arr == "syn-"],
    mapping = aes(x=ty, xend=ty, y=start, yend=end), inherit.aes = FALSE,
    arrow = arrow(15,length = unit(0.02,"npc"), type = "closed"),
    color = "red", size = 2
  ) +
  geom_segment(
    data = rbind(arr1,arr2)[arr == "ref"],
    mapping = aes(x=ty, xend=ty, y=start, yend=end), inherit.aes = FALSE,
    arrow = arrow(15,length = unit(0.02,"npc"), type = "closed"),
    color = "grey", size = 1
  ) +
#  geom_line(data = vac_only_reference) +
  geom_limits(limits[measure %in% c("combo.eff")]) +
  geom_path() +
  scale_size_vectorcontrol(guide=gds(1)) +
  scale_color_catchup(guide=gds(2)) +
  scale_linetype_vaccine(guide=gds(3)) +
  scale_alpha_manual("Calculated vs. Simulated", labels=c(obs="Full Simulation", est="Naive Estimation"), values = c(est=0.5, obs=1), guide=gds(4)) +
  scale_year() +
  scale_y_continuous(expand = expand_scale(0, 0)) +
  legend_theme

plotutil(p, h=5, w=7.5, args)