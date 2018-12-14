require(data.table)
require(ggplot2)

args <- c("foi_effstats.rds", "effstats.rds", "foi_fig_2.png")
args <- commandArgs(trailingOnly = TRUE)

base.stat.eff.dt <- readRDS(args[2])[vaccine == "cmdvi" & catchup == "none" & vc_coverage == 75]
base.stat.eff.dt[, foi := 1.0 ]

stat.eff.dt <- rbind(readRDS(args[1]), base.stat.eff.dt)



vec_lines <- c(`0`=3,`25`=2,`50`=5,`75`=1) # c(`25`="dotted",`50`="dashed",`75`="solid")

vartrans <- c(combo.eff="Combined Effectiveness", syn="Interaction")
varfact <- function(var) factor(vartrans[as.character(var)], levels = vartrans, ordered = T)

limits <- data.table(
  measure = rep("effectiveness", 4),
  vaccine = rep("none",4),
  vc_coverage = rep(0, 4),
  year = rep(1, 4),
  med = c(c(0, 1),c(-.15,.25))
)

no_vc <- stat.eff.dt[variable == "vac.eff",
  .(measure = "effectiveness", year, med=unique(med), vaccine, foi, vc_coverage = 0)
]

only_vc <- stat.eff.dt[variable == "vec.eff",
  .(measure = "effectiveness", year, med=unique(med), vaccine="none", foi, vc_coverage)
]

plot.dt <- rbind(stat.eff.dt[
  variable %in% c("combo.eff", "ind.eff"),
  .(measure = ifelse(variable=="combo.eff", "effectiveness", "pseudo"),
    year, med, vaccine, foi, vc_coverage)
], no_vc, only_vc)

p<-ggplot(plot.dt) +
  facet_grid(. ~ foi, scales="free_y", switch = "y") +
  aes(x=year, y=med, linetype=factor(vc_coverage), color=vaccine, alpha=factor(measure)) +
  geom_line() +
  geom_blank(data=limits) +
#  scale_color_manual("Vaccine Mech.", values=vac_cols) +
  scale_linetype_manual("VC Coverage %", values=vec_lines) +
  scale_x_continuous("Year", expand = expand_scale(0.03, 0)) +
  scale_y_continuous(expand = expand_scale(0, 0)) +
  theme_minimal() + theme(
    legend.direction = "horizontal",
    axis.title.y = element_blank(),
    strip.placement = "outside",
    legend.position = "bottom",
    legend.justification = c(0.5, 0.5),
    legend.margin = margin(), legend.spacing = unit(0, "pt"),
    panel.spacing.y = unit(15, "pt")
  ) +
  scale_alpha_manual(values=c("pseudo"=0.5,"effectiveness"=1))

ggsave(
  tail(args,1), p, device = "png",
  width = 7.5, height = 5, dpi = "retina", units = "in"
)