require(data.table)
require(lubridate)

args <- commandArgs(trailingOnly = TRUE)

## read in assorted input data
stat.eff10.dt <- readRDS(args[1])

# convenience function to look at sensitivity studies by appropriate
# dimensions -- uses "filt" to fix some params + pickout correct other layers
# first item is raw data (value = med.eff10, layer = "background")
# second item is smoothed data (value = smooth, layer = "foreground")
slice <- function(filt) rbind(stat.eff10.dt[
  eval(filt), .(doy, value=med.eff10, variable="Effectiveness",
    coverage, duration, durability,
    layer = "background"
  )], stat.eff10.dt[
  eval(filt), .(doy, value=smooth, variable="Effectiveness",
    coverage, duration, durability,
    layer = "foreground"
)])

# set duration, durability, floating coverage
coverage.dt <- slice(expression(duration == 90 & durability == 90))
# set coverage, durability, floating duration
duration.dt <- slice(expression(coverage == 75 & durability == 90))
# set duration, coverage, floating durability
durability.dt <- slice(expression(coverage == 75 & duration == 90))

saveRDS(coverage.dt, args[2])
saveRDS(duration.dt, args[3])
saveRDS(durability.dt, args[4])
