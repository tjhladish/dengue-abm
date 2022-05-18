library(DBI)
library(data.table)
library(lubridate)
library(ggplot2)
library(viridis)

setwd('~/documents/work/dengue/exp/tirs_study/')
db = dbConnect(RSQLite::SQLite(), "./tirs_trial-trans_loc_sensitivity.sqlite")
stmt = "select p.foi_mult, p.hh_mos_inc, m.* from par as p left join met as m on m.serial = p.serial;" # per_10k_AAR <-- metric for FOI
tirs_foi = dbGetQuery(db, stmt)
dbDisconnect(db)

tirs_foi = setDT(tirs_foi)
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