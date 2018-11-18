require(data.table)
require(ggplot2)

args <- c("effstats.rds","main_fig.png")
args <- commandArgs(trailingOnly = TRUE)

stat.eff.dt <- readRDS(args[1])

vac_mechs <- c("Traditional", "5th Sero")
vac_cols <- c("blue", "green")
names(vac_cols) <- vac_mechs
vacfact <- function(i) factor(vac_mechs[i+1], levels = vac_mechs, ordered = T)

vec_lines <- c(`0`=3,`25`=2,`50`=5,`75`=1) # c(`25`="dotted",`50`="dashed",`75`="solid")

catchuplabels <- c("No Catchup","Catchup")
cfact <- function(i) factor(catchuplabels[i+1], levels = catchuplabels, ordered = T)

vartrans <- c(combo.eff="Combined Effectiveness", syn="Interaction")
varfact <- function(var) factor(vartrans[as.character(var)], levels = vartrans, ordered = T)

limits <- data.table(
  variable = c(rep("combo.eff", 2), rep("syn", 2)),
  catchup = rep(0, 4),
  vac_mech = rep(0, 4),
  vc_coverage = rep(25, 4),
  year = rep(1, 4),
  med = c(c(0, 1),c(-.15,.25))
)

no_vc <- stat.eff.dt[variable == "vac.eff" & vc_coverage == 25,
  .(variable = "combo.eff", catchup, year, med, vac_mech, vc_coverage = 0)
]

p<-ggplot(stat.eff.dt[variable %in% c("combo.eff", "syn")]) +
  facet_grid(varfact(variable) ~ cfact(catchup), scales="free_y", switch = "y") +
  aes(x=year, y=med, color=vacfact(vac_mech), linetype=factor(vc_coverage)) +
  geom_line() +
  geom_line(data = no_vc) +
  geom_blank(data=limits) +
  scale_color_manual("Vaccine Mech.", values=vac_cols) +
  scale_linetype_manual("VC Coverage %", values=vec_lines) +
  theme_minimal() + theme(
    legend.direction = "horizontal",
    axis.title.y = element_blank(),
    strip.placement = "outside",
    legend.position = c(0.5, 0.5),
    legend.justification = c(0.5, 0.5),
    legend.margin = margin(), legend.spacing = unit(0, "pt")
  )

ggsave(
  tail(args,1), p, device = "png",
  width = 6, height = 6, dpi = "retina", units = "in"
)