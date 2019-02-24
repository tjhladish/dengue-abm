suppressPackageStartupMessages({
  require(data.table)
  require(cowplot)
})

args <- c("figref.rda", "rds/foi_effstats.rds", "rds/effstats.rds", "fig/fig_5.png")
args <- commandArgs(trailingOnly = TRUE)

tar <- tail(args, 1)

load(args[1])

base.stat.eff.dt <- readRDS(args[3])[vaccine == "t+cydtdv" & catchup == "routine" & vc_coverage == 75]
base.stat.eff.dt[, foi := 1.0 ]

stat.eff.dt <- rbind(
  readRDS(args[2]),
  base.stat.eff.dt
)[variable == "combo.eff"][, scenario := "vc+vac" ]

p<-ggplot(
  stat.eff.dt
) + theme_minimal() +
  aes(shape=vaccine, color=scenario, x=year+1, y=med, size=factor(vc_coverage)) +
  geom_ribbon(aes(color=NULL, fill=scenario, ymax=hi, ymin=lo), alpha=0.5) +
  geom_line() +
  geom_pchline(dt=stat.eff.dt, fill="white", show.legend = F) +
  scale_color_scenario(guide="none", aesthetics = c("color","fill")) +
  scale_size_vectorcontrol(guide="none") +
  scale_shape_vaccine(guide="none") +
  scale_year() +
  scale_effectiveness() +
  facet_grid(. ~ foi, labeller=facet_labels) +
  FOIfacettitle +
  coord_cartesian(clip="off", ylim=c(-.125,1), xlim=c(0,40)) + theme(
    panel.spacing.x = unit(12, "pt")
  )

plotutil(p, h=3, w=6.75, tar)