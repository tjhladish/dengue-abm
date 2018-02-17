require(reshape2)
require(data.table)
require(lubridate)

args <- commandArgs(trailingOnly = TRUE)

weather.dt = fread(args[1],
 select = c("NAME","DATE","PRCP"),
 colClasses = c(NAME="factor", PRCP="integer"))[ # per ?fread: "Dates are read as character currently."
   NAME=='AEROP.INTERNACIONAL'
][, DATE := as.Date(DATE) ][, rained := PRCP > 0 ]

rain.dt <- weather.dt[,
  .(prob = mean(rained, na.rm=T)), keyby=.(month=month(DATE), day=day(DATE))
][ !(month == 2 & day == 29) ][, doy := .I ]
# drop leap year extra day


# smooth data similar to R0, but we need the values shifted earlier by one week (7 days)
# lag by one week to account for egg-to-emergence time, based on Table 2 in Rueda et al, Temperature-dependent
# development and survival rates of Culex quinquefasciatus and Aedes aegypti (Diptera: Culicidae), J Med Ent, 1990, 27:5
rain.dt[, smooth := smooth.spline(rep(prob,3), spar=0.6)$y[1:.N+.N-7] ]

store <- rbind(
  rain.dt[,.(doy, value = smooth / max(smooth), variable = "Mos. pop.",
    coverage = "reference", duration = "reference", durability = "reference",
    layer = "foreground"
  )],
  rain.dt[,.(doy, value = prob / max(smooth), variable = "Mos. pop.",
    coverage = "reference", duration = "reference", durability = "reference",
    layer = "background"
  )]
)

saveRDS(store, args[2])
