suppressPackageStartupMessages({
  require(data.table)
  require(cowplot)
})

args <- c("figref.rda", "rds/foi_effstats.rds", "rds/effstats.rds", "fig/fig_5.png")
args <- commandArgs(trailingOnly = TRUE)

tar <- tail(args, 1)

load(args[1])

base.stat.eff.dt <- readRDS(args[3])[vaccine == "cmdvi" & catchup == "routine" & vc_coverage == 75]
base.stat.eff.dt[, foi := 1.0 ]

stat.eff.dt <- rbind(readRDS(args[2]), base.stat.eff.dt)[variable %in% c("combo.eff","vac.eff","vec.eff","ind.eff")]
stat.eff.dt[,
  scenario := ifelse(variable == "combo.eff", "vc+vac",
              ifelse(variable == "vac.eff"  , "vac",
              ifelse(variable == "vec.eff"  , "vc", "not")))
]
stat.eff.dt[scenario == "vac", vc_coverage := 0]
stat.eff.dt[scenario == "vc", vaccine := "none"]
# vaccine line
# vc line
# combo line

real.dt <- stat.eff.dt[variable %in% c("combo.eff","vac.eff","vec.eff")]

# annos <- list(
#   annotate("line", x=naive.eff.cmdvi$year+1, y=naive.eff.cmdvi$value, size=vc_sizes["75"], linejoin = "mitre", lineend = "butt", color = "#AAAAFF"),
#   annotate("point", x=naive.eff.cmdvi.thin$year+1, y=naive.eff.cmdvi.thin$value, shape=vac_pchs["cmdvi"], size=2, color="#55CC55", fill="#55CC55")
# )

p<-ggplot(
  real.dt
) + theme_minimal() +
  aes(shape=vaccine, color=scenario, x=year+1, y=med, size=factor(vc_coverage)) +
  geom_line(data=stat.eff.dt[variable == "ind.eff"], color="#AAAAFF") +
  geom_point(data=stat.eff.dt[variable == "ind.eff"][((year+1) %% 5 == 0) | year == 0], size=2, fill="white", color="#55CC55") +
  geom_line() +
  geom_point(data=real.dt[((year+1) %% 5 == 0) | year == 0], size=2, fill="white") +
  scale_color_scenario(guide="none") +
  scale_size_vectorcontrol(guide="none") +
  scale_shape_vaccine(guide="none") +
  scale_year() +
  scale_effectiveness() +
  facet_grid_freey(. ~ foi, labeller=facet_labels) +
  coord_cartesian(clip="off", ylim=c(-.125,1), xlim=c(0,40)) + theme(
    panel.spacing.x = unit(12, "pt")
  )

plotutil(p, h=3, w=6.75, tar)