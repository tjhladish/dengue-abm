require(data.table)
require(RSQLite)
require(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
# args <- c("~/Dropbox/who/irs_results-refit2-partial.sqlite")

drv = dbDriver("SQLite")

db = dbConnect(drv, args[1], flags=SQLITE_RO)

plan.dt <- data.table(dbGetQuery(db,
  "SELECT vector_control,
   timing+1 AS doy, vc_coverage*100 AS coverage,
   campaign_duration AS duration, eff_days AS durability,
   status,
   posterior AS particle
   FROM par P
   JOIN job J ON J.serial = P.serial
   WHERE strat_years != 10;" # ignore the stopping scenarios
))

dbDisconnect(db)

plan.dt[vector_control==0, coverage := 0]
plan.dt[, cov := sprintf("%i%% hh cov.", coverage)]
plan.dt[, scenario := ifelse(vector_control, "intervention","baseline")]
plan.dt[, rollout := factor(ifelse(duration == 0, "1 day roll", ifelse(duration == 1, "90 day roll", "continuous")))]
plan.dt[, dur := sprintf("%i day eff.", durability)]

kb <- grep("particle", names(plan.dt), invert=T, value=T)

res = 300
mag = 0.85*res/72
png(args[2], width = 1000*mag, height = 1000*mag, units = "px", res=res)
ggplot(plan.dt[, .N, keyby=kb]) + theme_minimal() +
  aes(x = doy, y = N, fill=status) +
  facet_grid(scenario + cov ~ rollout + dur) +
  geom_col() +
  scale_fill_manual(values=c(D="green",Q="red",R="blue", P="black")) +
  scale_x_continuous(name="start day-of-year") +
  scale_y_continuous(name="# of simulation runs")
dev.off()