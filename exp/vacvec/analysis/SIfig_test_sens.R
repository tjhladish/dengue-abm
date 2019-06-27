suppressPackageStartupMessages({
  require(data.table)
  require(cowplot)
})

args <- c("figref.rda", "rds/nolag_effstats.rds", "fig/SIfig_9.png")
args <- commandArgs(trailingOnly = TRUE)

tar <- tail(args, 1)

load(args[1])
effstats.dt <- readRDS(args[2])[
  vaccine != "d70e" & foi == 1.0 & variable %in% c("combo.eff","c.combo.eff","vac.eff","c.vac.eff") & vc_coverage == 75
]

plotter <- function(dt, efflabel) ggplot(dt[!is.na(false_neg) & false_neg != 0]) + 
  aes(year, med, linetype=catchup) +
  facet_grid(false_neg ~ false_pos, labeller = labeller(
    false_neg = function(p) sprintf("FNR=%s", p),
    false_pos = function(p) sprintf("FPR=%s", p)
  )) +
  geom_line(aes(color="none", alpha="background"), dt[is.na(false_neg), .(year, med, catchup)]) +
  geom_line(aes(color="perfect", alpha="background"), dt[false_neg==0, .(year, med, catchup)]) +
  geom_line(aes(color="realistic", alpha="foreground")) + scale_year() + scale_effectiveness(efflabel) +
  scale_color_manual("Testing", values=c(none="firebrick",perfect="dodgerblue",realistic="black")) +
  scale_alpha_manual(values=c(background=0.5, foreground=1), guide="none") +
  theme_minimal()

annual <- effstats.dt[variable == "combo.eff"]

annualp <- plotter(annual, "Annual Combined Effectiveness")

cumul <- effstats.dt[variable == "c.combo.eff"]

cumulp <- plotter(cumul, "Cumulative Combined Effectiveness")

annual.vac <- effstats.dt[variable == "vac.eff"]
pannualv <- plotter(annual.vac, "Vaccine-Only Annual Effectiveness")
cumul.vac <- effstats.dt[variable == "c.vac.eff"]
pcumulv <- plotter(cumul.vac, "Vaccine-Only Cumulative Effectiveness")
