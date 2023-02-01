library(DBI)
library(data.table)
library(lubridate)
library(ggplot2)
library(viridis)

if (interactive()) { setwd("~/documents/work/dengue/exp/tirs_study") }

.args <- if (interactive()) c(
  "tirs_trial-trans_loc_sensitivity-v2.0.sqlite"
) else commandArgs(trailingOnly = TRUE)

db_path <- .args[1]

exp_str <- sub("(*).sqlite", "\\1", db_path)
figPath <- paste0("./fig/", exp_str, "_foi_plots")
dir.create(figPath, recursive = T)

db = dbConnect(RSQLite::SQLite(), db_path)
stmt = "select p.foi_mult, p.hh_mos_inc, m.* from par as p left join met as m on m.serial = p.serial;" # per_10k_AAR <-- metric for FOI
tirs_foi = dbGetQuery(db, stmt)
dbDisconnect(db)

tirs_foi = setDT(tirs_foi)

tirs_foi[, per_10k_AAR_normd := per_10k_AAR/mean(per_10k_AAR[foi_mult==1.0]), by = .(hh_mos_inc)]
tirs_foi[, `:=`(
  eff_total_inf = 1-((inf_arm2_y1 + inf_arm2_y2)/(inf_arm1_y1 + inf_arm1_y2)),
  eff_total_case = 1-((sym_arm2_y1 + sym_arm2_y2)/(sym_arm1_y1 + sym_arm1_y2)),
  eff_total_dss = 1-((dss_arm2_y1 + dss_arm2_y2)/(dss_arm1_y1 + dss_arm1_y2))
)]

inf.quant <- tirs_foi[, .(
  lower=quantile(eff_total_inf, c(0.25), na.rm=T, names=F), 
  median=quantile(eff_total_inf, c(0.5), na.rm=T, names=F),
  upper=quantile(eff_total_inf, c(0.75), na.rm=T, names=F)
  ), by = .(foi_mult, hh_mos_inc)]

case.quant <- tirs_foi[, .(
  lower=quantile(eff_total_case, c(0.25), na.rm=T, names=F), 
  median=quantile(eff_total_case, c(0.5), na.rm=T, names=F),
  upper=quantile(eff_total_case, c(0.75), na.rm=T, names=F)
), by = .(foi_mult, hh_mos_inc)]

dss.quant <- tirs_foi[, .(
  lower=quantile(eff_total_dss, c(0.25), na.rm=T, names=F), 
  median=quantile(eff_total_dss, c(0.5), na.rm=T, names=F),
  upper=quantile(eff_total_dss, c(0.75), na.rm=T, names=F)
), by = .(foi_mult, hh_mos_inc)]

# scatter plot of eff vs. per_10k_AAR
inf.plot <- ggplot(inf.quant) +
  geom_line(aes(foi_mult, median)) +
  geom_point(aes(foi_mult, median)) +
  geom_ribbon(aes(foi_mult, ymin = lower, ymax = upper), alpha = 0.3) +
  facet_wrap(vars(hh_mos_inc)) +
  theme_light() +
  labs(x = "FOI multiplier", y = "TIRS effectiveness vs. infections")

case.plot <- ggplot(case.quant) +
  geom_line(aes(foi_mult, median)) +
  geom_point(aes(foi_mult, median)) +
  geom_ribbon(aes(foi_mult, ymin = lower, ymax = upper), alpha = 0.3) +
  facet_wrap(vars(hh_mos_inc)) +
  theme_light() +
  labs(x = "FOI multiplier", y = "TIRS effectiveness vs. cases")

dss.plot <- ggplot(dss.quant) +
  geom_line(aes(foi_mult, median)) +
  geom_point(aes(foi_mult, median)) +
  geom_ribbon(aes(foi_mult, ymin = lower, ymax = upper), alpha = 0.3) +
  facet_wrap(vars(hh_mos_inc)) +
  theme_light() +
  labs(x = "FOI multiplier", y = "TIRS effectiveness vs. DSS")

aar_v_foi.plot <- ggplot(tirs_foi) +
  geom_violin(aes(x = factor(foi_mult), y = per_10k_AAR, fill = factor(foi_mult))) +
  facet_wrap(vars(hh_mos_inc)) +
  scale_fill_viridis_d() +
  theme_light() +
  theme(legend.position = "none") +
  labs(x = "FOI multiplier", y = "Annual attack rate per 10k")

ggplot(tirs_foi[is.finite(eff_total_inf)]) +
  geom_point(aes(x = per_10k_AAR_normd, y = eff_total_inf), alpha = 0.1) +
  facet_wrap(vars(hh_mos_inc)) +
  ylim(c(0, 1))

ggsave(filename = paste0(figPath, "/", exp_str, "_inf_eff.png"), plot = inf.plot, device = 'png', units = 'in', height = 8, width = 10, dpi = 300)
ggsave(filename = paste0(figPath, "/", exp_str, "_case_eff.png"), plot = case.plot, device = 'png', units = 'in', height = 8, width = 10, dpi = 300)
ggsave(filename = paste0(figPath, "/", exp_str, "_dss_eff.png"), plot = dss.plot, device = 'png', units = 'in', height = 8, width = 10, dpi = 300)
ggsave(filename = paste0(figPath, "/", exp_str, "_aar_v_foi.png"), plot = aar_v_foi.plot, device = 'png', units = 'in', height = 8, width = 10, dpi = 300)
