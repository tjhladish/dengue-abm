require(data.table)
require(ggplot2)

source("projref.R")

args <- c("effstats.rds","fig_2.png")
args <- commandArgs(trailingOnly = TRUE)

stat.eff.dt <- readRDS(args[1])

vec_lines <- c(`0`=3,`25`=2,`50`=5,`75`=1) # c(`25`="dotted",`50`="dashed",`75`="solid")

vartrans <- c(combo.eff="Combined Effectiveness", syn="Interaction")
varfact <- function(var) factor(vartrans[as.character(var)], levels = vartrans, ordered = T)

limits <- data.table(
  variable = c(rep("combo.eff", 2), rep("syn", 2)),
  catchup = rep(reference.scenario$catchup[1], 4),
  vaccine = rep(reference.scenario$vaccine[1], 4),
  vc_coverage = rep(25, 4),
  year = rep(1, 4),
  med = c(c(0, 1),c(-.15,.25))
)

no_vc <- stat.eff.dt[variable == "vac.eff",
  .(variable = "combo.eff", catchup, year,
    med={
      if (length(unique(med)) != 1) warning("non unique median vac.eff")
      median(med)
    },
    vaccine, vc_coverage = 0)
]

p<-ggplot(stat.eff.dt[variable %in% c("combo.eff", "syn")]) +
  facet_grid(varfact(variable) ~ catchup, scales="free_y", switch = "y") +
  aes(x=year, y=med, color=vaccine, linetype=factor(vc_coverage)) +
  geom_line() +
  geom_line(data = no_vc) +
  geom_blank(data=limits) +
  scale_color_manual("Vaccine Mech.", values=vac_cols) +
  scale_linetype_manual("VC Coverage %", values=vec_lines) +
  scale_x_continuous("Year", expand = expand_scale(0.03, 0)) +
  scale_y_continuous(expand = expand_scale(0, 0)) +
  theme_minimal() + theme(
    legend.direction = "horizontal",
    axis.title.y = element_blank(),
    strip.placement = "outside",
    legend.position = c(0.5, 0.5),
    legend.justification = c(0.5, 0.5),
    legend.margin = margin(), legend.spacing = unit(0, "pt"),
    panel.spacing.y = unit(15, "pt")
  )

ggsave(
  tail(args,1), p, device = "png",
  width = 7.5, height = 5, dpi = "retina", units = "in"
)