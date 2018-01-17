require(data.table)

args <- commandArgs(trailingOnly = TRUE)

sim.raw = fread(args[1], header=T, col.names=c('serial','day','total_cases'))
sim.raw = sim.raw[day<=50004,] # end before 2016
# #sim.raw = sim.raw[day>46354,]  # start with 2006

#sim.raw[, inc := intro + local ]
sim.raw[, doy := day %% 365 + 1 ]

store <- sim.raw[, .(
  value = mean(total_cases),#value = mean(inc),
  variable = "Cases (model)",
  coverage = "reference",
  duration = "reference",
  durability = "reference",
  layer = "foreground"
), by=doy]

store[, value := value/max(value) ]

saveRDS(store, args[2])
