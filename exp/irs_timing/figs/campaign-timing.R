rm(list=ls())

args <- commandArgs(trailingOnly = T)
# args <- paste0("~/Dropbox/", c("who/irs_timing-summer_winter-foi_refit2.sqlite", "daily-irs_refit2-simple-uniq.out", "who/fig1_data/campaign-timing-prevalence.local.rds"))
require(RSQLite)
require(data.table)

drv = dbDriver("SQLite")
#db.path <- sprintf("%s%s/%s",.dbox.path,"who","irs_timing-summer_winter-foi.sqlite")
#log.path <- sprintf("%s%s",.dbox.path,"daily-irs_refit-simple-nofilenames-uniq.out")

db = dbConnect(drv, args[1])

# campaign_duration is always 1
# coverage is always 75%

foilvls <- c("70% Mosquito Population", "Baseline Mosquito Population", "130% Mosquito Population")
timelvls <- c("None", "Proactive", "Reactive")

ref.dt <- data.table(
  dbGetQuery(db, sprintf(
    'SELECT P.serial AS serial,
       CASE WHEN vector_control = 0 THEN "%s" WHEN timing = 147 THEN "%s" ELSE "%s" END AS intervention,
       CASE WHEN foi < 1.0          THEN "%s" WHEN foi = 1.0    THEN "%s" ELSE "%s" END AS foi
     FROM par P, job J
       WHERE P.serial = J.serial
         AND status = \'D\';',
      timelvls[1], timelvls[2], timelvls[3],
      foilvls[1], foilvls[2], foilvls[3]
  )), stringsAsFactors = T
)

dbDisconnect(db)

ref.dt[,
  foi := factor(foi, levels = foilvls, ordered = T)
][,
  intervention := factor(intervention, levels = timelvls, ordered = T)  
]

target <- gsub(".+-(\\w+\\.\\w+).rds$","\\1",args[3])
tarcols <- c("incidence.intro","incidence.local","prevalence.intro","prevalence.local")
trans <- sprintf("V%d", 1:length(tarcols) + 2)
dropcols <- trans[which(tarcols != target)]
## TODO check for .Rdata

daily.dt <- fread(
  args[2], sep = " ", header = F,
  col.names = c("serial","day", "value"),
  drop = dropcols
)[ref.dt, on="serial"][, day    := day - min(day) ]

plot.dt <- daily.dt[,{
  mn <- mean(value)
  ps <- quantile(value, probs = c(0,.25,.5,.75, 1))
  .(ave=mn, min=ps[1], lo=ps[2], md=ps[3], hi=ps[4], max=ps[5])
  }, keyby=.(foi, intervention, day)
]

saveRDS(plot.dt, args[3])