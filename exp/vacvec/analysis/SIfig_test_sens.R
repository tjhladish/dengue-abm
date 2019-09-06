suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
  require(cowplot)
})

args <- c("figref.rda", "rds/nolag_effstats.rds", "fig/SIfig_9a.png")
args <- commandArgs(trailingOnly = TRUE)

tar <- tail(args, 1)

load(args[1])
effstats.dt <- readRDS(args[2])[
  vaccine != "d70e" & foi == 1.0 & variable %in% c("combo.eff","c.combo.eff","vac.eff","c.vac.eff") & vc_coverage == 75
]

plotter <- function(dt, efflabel) ggplot(dt[!is.na(false_neg) & false_neg != 0]) + 
  aes(year, med, linetype=catchup) +
  facet_grid(false_neg ~ false_pos, labeller = labeller(
    false_neg = function(p) sprintf("%i%% False Sero-", as.numeric(p)*100),
    false_pos = function(p) sprintf("%i%% False Sero+", as.numeric(p)*100)
  )) +
  geom_line(aes(color="none", alpha="background"), dt[is.na(false_neg), .(year, med, catchup)]) +
  geom_line(aes(color="perfect", alpha="background"), dt[false_neg==0, .(year, med, catchup)]) +
  geom_line(aes(color="realistic", alpha="foreground")) + scale_year() + scale_effectiveness(efflabel) +
  coord_cartesian(ylim=c(0,1)) +
  scale_color_manual("Testing",
    labels = c(none="None",perfect="Perfect",realistic="Imperfect"),
    values=c(none="firebrick",perfect="dodgerblue",realistic="black"),
    guide = guide_legend(title.position = "top")
  ) +
  scale_alpha_manual(values=c(background=0.5, foreground=1), guide="none") +
  scale_linetype_manual("Vaccine Program",
    labels = c(`vc+vac`="Routine w/ Catchup", `routine`="Routine-Only"),
    values = c(`vc+vac`="solid",`routine`="dashed"),
    guide = guide_legend(title.position = "top")
  ) +
  theme_minimal() + theme(
    legend.margin = margin(l=10, r=10),
    panel.spacing.y = unit(2,"line"),
    legend.title = element_text(size=rel(0.6), hjust = 0.5),
    legend.text = element_text(size=rel(0.6)),
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.key.height = unit(10,"pt"),
    legend.position = c(0.5, 0.5)
  )

annual <- effstats.dt[variable == "combo.eff"]

annualp <- plotter(annual, "75% TIRS & CYD-TDV Annual Effectiveness")

saveplotter <- function(p, sb) save_plot(gsub("a",sb,tar), p, ncol = 2, nrow = 2, base_height = 2.5, base_aspect_ratio = 1.1)

saveplotter(annualp, "a")

cumul <- effstats.dt[variable == "c.combo.eff"]

cumulp <- plotter(cumul, "75% TIRS & CYD-TDV Cumulative Effectiveness")

saveplotter(cumulp, "b")

annual.vac <- effstats.dt[variable == "vac.eff"]
pannualv <- plotter(annual.vac, "CTD-TDV Only Annual Effectiveness")

saveplotter(pannualv, "c")

cumul.vac <- effstats.dt[variable == "c.vac.eff"]
pcumulv <- plotter(cumul.vac, "CYD-TDV Only Cumulative Effectiveness")

saveplotter(pcumulv, "d")