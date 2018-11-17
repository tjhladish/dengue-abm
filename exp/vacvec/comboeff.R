require(data.table)

args <- c("effectiveness.rds", "comboeff.rds")
args <- commandArgs(trailingOnly = TRUE)
# args <- c("~/Dropbox/irs_timing-refit0_intro-fix.sqlite", "~/Dropbox/who/fig1_data/baseline.rds")

effectiveness.dt <- readRDS(args[1])

vec.dt <- effectiveness.dt[vac == 0, .(vec.eff=eff), by=.(particle, vc_coverage, year)]
vac.dt <- effectiveness.dt[vc == 0, .(vac.eff=eff), by=.(particle, vac_mech, catchup, year)]
combo.dt <- effectiveness.dt[(vc != 0) & (vac != 0)]

syn.dt <- combo.dt[
  vec.dt, on=.(particle, vc_coverage, year) # join vector control only results
][
  vac.dt, on=.(particle, vac_mech, catchup, year) # join vaccine only results
][, # get the interesting measures
  .(
    combo.eff=eff, ind.eff = (vec.eff + vac.eff - vec.eff*vac.eff),
    vec.eff, vac.eff
  ), by=.(particle, year, vc_coverage, vac_mech, catchup)
  # ...organized by relevant divisions
]

syn.dt[,
  syn := combo.eff - ind.eff
][,
  # if syn positive, should be fraction of space between 1 and ind.eff
  syn.frac := syn / ifelse(syn < 0, ind.eff, 1-ind.eff)
]

saveRDS(syn.dt, args[2])

stat.syn.dt <- syn.dt[,
  .(
    med.syn = stats::median(syn, na.rm=T),
    med.syn.frac = stats::median(syn.frac, na.rm=T)
  ),
  #.(med.eff = sum(eff)),
  keyby=.(
    vc_coverage,
    vac_mech,
    catchup,
    year
  )
]

require(ggplot2)

vac_mechs <- c("5th Sero", "Traditional")
vac_cols <- c("green","blue")
names(vac_cols) <- vac_mechs

vec_lines <- c(`25`="dotted",`50`="dashed",`75`="solid")

ggplot(stat.syn.dt) + aes(
  y = med.syn.frac, x = year,
  color = factor(vac_mechs[vac_mech+1]),
  linetype = factor(vc_coverage)
) + facet_grid(. ~ factor(c("No Catchup","Catchup")[catchup+1])) +
  geom_line() +
  scale_color_manual("Vaccine Mech.", values=vac_cols) +
  scale_linetype_manual("VC Coverage %", values=vec_lines) +
  theme_minimal() + theme(
    legend.direction = "horizontal",
    legend.position = c(0.5, 1),
    legend.justification = c(0, 1),
    legend.margin = margin(), legend.spacing = unit(0, "pt")
  ) + labs(y="Interaction Fraction ([-1,1] scale)")

ggplot(stat.syn.dt) + aes(
  y = med.syn, x = year,
  color = factor(vac_mechs[vac_mech+1]),
  linetype = factor(vc_coverage)
) + facet_grid(. ~ factor(c("No Catchup","Catchup")[catchup+1])) +
  geom_line() +
  scale_color_manual("Vaccine Mech.", values=vac_cols) +
  scale_linetype_manual("VC Coverage %", values=vec_lines) +
  theme_minimal() + theme(
    legend.direction = "horizontal",
    legend.position = c(0.5, 1),
    legend.justification = c(0, 1),
    legend.margin = margin(), legend.spacing = unit(0, "pt")
  ) + labs(y="Absolute Interaction Delta")