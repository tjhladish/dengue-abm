suppressPackageStartupMessages({
  require(data.table)
  require(cowplot)
}) 

warnnonunique <- function(var, variable, collapse = median) {
  if (length(unique(var)) != 1) warning(sprintf("non unique %s", unique(variable)))
  collapse(var)
}

args <- c("figref.rda", "rds/effstats.rds","fig/fig_4.png")
args <- commandArgs(trailingOnly = TRUE)

load(args[1])
effstats.dt <- readRDS(args[2])
tar <- tail(args, 1)

cmb.eff <- effstats.dt[variable %in% c("combo.eff","c.combo.eff"), .(
  value = warnnonunique(med, variable),
  lo=warnnonunique(lo, variable),
  hi=warnnonunique(hi, variable),
  scenario = trans_scnario(1, 1),
  estimate = "simulated"
), keyby=.(
  vaccine = factor(vaccine, rev(levels(vaccine)), ordered = T), catchup, vc_coverage, year,
  measure = factor(gsub("combo.","",variable,fixed = T), levels=c("eff","c.eff"), ordered = T)
)]

p <- ggplot(cmb.eff) + aes(
  shape=vaccine, color=scenario, size=factor(vc_coverage),
  x=year+1, y=value, group=interaction(vaccine, vc_coverage, catchup, measure)
) + theme_minimal() +
  facet_grid_freey(measure ~ vc_coverage, labeller = facet_labels) +
  geom_ribbon(
    aes(fill=scenario, color=NULL, ymin=lo, ymax=hi),
    alpha=0.5, show.legend = F
  ) +
  geom_line(alpha=1, size=vc_sizes["0"]) +
  geom_pchline(dt=cmb.eff[catchup == "routine"], fill="white", alpha=1) +
  geom_pchline(dt=cmb.eff[catchup != "routine"], fill="black", alpha=1) +
  scale_year() + scale_effectiveness() +
  # scale_fill_interaction(
  #   guide = gds(1, keyheight=unit(12,"pt"), label.position = "right", direction="vertical", override.aes=list(alpha=c(0.4,0.4)))
  # ) +
  scale_pchlty_vaccine(guide = "none") +
  scale_color_scenario(guide = "none", value="black", aesthetics = c("color","fill")) +
  scale_size_vectorcontrol(guide="none") +
  coord_cartesian(ylim=c(0,1), xlim=c(0,40), clip="off") +
  TIRSfacettitle +
  theme(
    legend.margin = margin(), legend.spacing = unit(25, "pt"),
    legend.text = element_text(size=rel(0.5)),
    legend.title = element_text(size=rel(0.6)), legend.title.align = 0.5,
    panel.spacing.y = unit(15, "pt"), panel.spacing.x = unit(15, "pt"),
    legend.key.height = unit(1,"pt"),
    legend.box.spacing = unit(2.5, "pt"),
    legend.position = c(100/120, 0.6) # think this looks best, but can comment out to return to margin
  ) +
  scale_alpha_manual(values=c(delta=int_alpha), guide = "none")

plotutil(p, h=5, w=7.5, tar)