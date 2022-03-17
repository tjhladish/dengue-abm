
require(data.table)

## get script args; for debugging, uncomment section that follows
.args <- if (interactive()) c(
  file.path("~", "Dropbox", "who", "fig1_data",
    c(
      sprintf("%s.rds", c("obs", "sim", "eip", "R0", "mos")),
      "fig1.csv"
    )
  )
) else commandArgs(trailingOnly = TRUE)

## read in assorted input data
obs.dt <- readRDS(.args[1])[, .(doy, value, variable = "mean_obs_cases_frac_max")]
sim.dt <- readRDS(.args[2])[, .(doy, value, variable = "mean_sim_cases_frac_max")]
eip.dt <- readRDS(.args[3])[, .(doy, value, variable = "eip_rel_max")]
R0.dt  <- readRDS(.args[4])[, .(doy, value, variable = "transmission_multiplier")]
mos.dt <- readRDS(.args[5])[, .(doy, value, variable = "mos_pop_rel_max")]

fwrite(rbind(
  obs.dt, sim.dt, eip.dt, R0.dt, mos.dt
), tail(.args, 1))
