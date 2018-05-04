require(data.table)
require(lubridate)

args <- commandArgs(trailingOnly = TRUE)
## args <- c("~/Dropbox/who/fig1_data/stat.eff10.rds")

## read in assorted input data
stat.eff10.dt <- readRDS(args[1])

# convenience function to look at sensitivity studies by appropriate
# dimensions -- uses "filt" to fix some params + pickout correct other layers
# first item is raw data (value = med.eff10, layer = "background")
# second item is smoothed data (value = smooth, layer = "foreground")
slice <- function(filt) rbind(stat.eff10.dt[
  eval(filt), .(doy, value=med.eff10, variable="Effectiveness",
    coverage, duration, durability, efficacy,
    layer = "background"
  )], stat.eff10.dt[
  eval(filt), .(doy, value=smooth, variable="Effectiveness",
    coverage, duration, durability, efficacy,
    layer = "foreground"
)])

## TODO: introspect sensitivities
#   have sensitivity names;
#   can get (1) the main value for those

# set duration, durability, efficacy, floating coverage
coverage.dt <- slice(expression(duration == 90 & durability == 90 & efficacy == 0.8))
# set coverage, durability, efficacy, floating duration
duration.dt <- slice(expression(coverage == 75 & durability == 90 & efficacy == 0.8))
# set duration, coverage, efficacy, floating durability
durability.dt <- slice(expression(coverage == 75 & duration == 90 & efficacy == 0.8))
# set duration, coverage, durability, floating efficacy
efficacy.dt <- slice(expression(coverage == 75 & duration == 90 & durability == 90))

saveRDS(coverage.dt, args[2])
saveRDS(duration.dt, args[3])
saveRDS(durability.dt, args[4])
saveRDS(efficacy.dt, args[5])
