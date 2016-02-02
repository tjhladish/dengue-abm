#!/usr/bin/env Rscript
# must be run in aggregated files directory
tar <- "~/Downloads/who-feb-2016/who-feb-2016-aggregated/"
setwd(tar)

rm(list=ls())

args <- commandArgs(trailingOnly = T)
poppath <- ifelse(is.na(args[1]), "~/git/dengue/pop-merida/pop-merida/population-merida.txt", args[1])
dbpath <- ifelse(is.na(args[2]), "~/Dropbox", args[2])

cat("loading pop from: ",poppath,"\n")
cat("dropbox path\n", dbpath,"\n")

srcfiles <- list.files(pattern = "^[0-4].+\\..+")

if(length(srcfiles) == 0) stop("no source files")

require(reshape2)
require(data.table)
require(parallel)

#  scenario [0-4][01]{6}
#  1transmission level
#  2catchup?
#  3vac mechanism
#  4vaccination?
#  5target age
#  6catch age
#  7coverage

translate_scenario <- function(dt) { # assumes transmission already separated
  ky <- key(dt)
  dt[substr(scenario,6,6) == "0", scen := "coverageAT50%"]
  dt[substr(scenario,2,2) == "1", scen := "altVaccine" ]
  dt[substr(scenario,1,1) == "1", scen := "catchUp"]
  dt[substr(scenario,5,5) == "1", scen := paste0(scen,"To30")]
  dt[grepl("^0{2}10{2}1",scenario), scen := "reference"]
  dt[substr(scenario,4,4) == "1", scen := "routineAT16"]
  dt[substr(scenario,3,3) == "0", scen := "noVaccine"]
  setkeyv(subset(dt[, scenario:=scen ], select=-scen), ky)
}

rr_severity_sec_vs_pri = 20
ph1 = 0.111
ph2 = ph1*1.88
phs = 1
# ps2 = rr_severity_sec_vs_pri*ps1
# ph1 = phs*ps1 + phc*pc1
# ph2 = phs*ps2 + phc*pc2
# 
# 0.111 = phs*ps1 + phc*(1-ps1)
# 0.111*1.88 = phs*20*ps1 + phc*(1-20*ps1)
# 0.111*1.88 = 20*ps1 + phc*(1-20*ps1)
# 0.111*(1-1.88) = -19*ps1 - 0.111*1.88*ps1 + 20*0.111*ps1
# 0.111*(1-1.88) = (-19 - 0.111*1.88 + 20*0.111)*ps1
ps1 = ph1*(1-1.88) / (-(rr_severity_sec_vs_pri-1) - ph1*1.88 + rr_severity_sec_vs_pri*ph1)
# ps2 = 20*ps1 = 0.1149942
# pc1 = 1 - ps1
# pc2 = 1 - ps2
phc = (0.111 - ps1)/(1-ps1)

# require(ggplot2)

cfr <- .00078

age_cats <- c("<9yrs","9-18yrs","19+yrs","overall")

pop <- fread(poppath)[,list(pop=.N),keyby=age][,
  age_category := factor(ifelse(
    age < 9, age_cats[1], ifelse(
      age < 19, age_cats[2],
      age_cats[3]
    )), levels = age_cats, ordered = TRUE)
]

src <- subset(merge(
  dcast.data.table(melt(
    rbindlist(lapply(srcfiles, fread, colClasses = c(scenario="character"))),
    id.vars = c("serial","scenario","age","outcome"), variable.name = "year", value.name = "count"
  )[,
    year := as.integer(gsub("y", "", year))
  ][,
    outcome := factor(c("infection","mild","severe")[outcome+1], levels=c("infection","mild","severe"), ordered = T)
  ], serial + scenario + age + year ~ outcome, value.var = "count"
), pop, by = "age"), select=-age_category)[,
    particle_id := floor(serial/80)
  ][,
    `:=`(
      infections = infection + mild + severe,
      symptomatic_cases = mild + severe,
      hospitalised_cases = mild*phc + severe
    ) 
  ][,
    deaths := hospitalised_cases * cfr
  ][,
    `:=`(
      transmission_setting = seq(10,90,20)[as.integer(substr(scenario, 1, 1))+1],
      scenario = substr(scenario,2,7)
    )
  ]

translate_scenario(src)
slc <- function(dt, ...) dt[...,
  list(infections, symptomatic_cases, hospitalised_cases, deaths),
  keyby=list(age, year, pop, particle_id, transmission_setting, scenario)
]

ref <- subset(slc(src, scenario == "noVaccine"), select=-scenario)
forecasts <- slc(src, scenario != "noVaccine")

mvars = c("infections", "symptomatic_cases", "hospitalised_cases", "deaths")
mlt.ref <- melt(ref, measure.vars = mvars, variable.name = "outcome")
setkeyv(mlt.ref, grep("value", names(mlt.ref), value = T, invert = T))
mlt.fore <- setkeyv(
  melt(forecasts, measure.vars = mvars, variable.name = "outcome"),
  c(key(mlt.ref),"scenario")
)

agetrans <- pop[,age_category,keyby=age]
cached_res <- mlt.fore[mlt.ref][agetrans] # cached_res <- readRDS("~/Downloads/who-rawpre.rds")
setkeyv(cached_res, c(key(cached_res), "age_category"))

cumulative_events <- setkeyv(cached_res[,
  list(year=paste0("cum", year+1), value=cumsum(value), i.value=cumsum(i.value)),
  keyby=eval(grep("year", key(cached_res), invert = T, value=T))
][, averted := i.value - value ][, outcome := paste0(outcome, "_averted")], key(cached_res))

tarevents <- cumulative_events[year %in% c("cum10","cum30")][,
  list(div=sum(i.value), averted=sum(averted), pop=sum(pop)),
  keyby=eval(c(grep("(age|pop)",key(cumulative_events),invert=T,value=T),"age_category"))
]

tareventsoverall <- tarevents[,
  list(averted=sum(averted), div=sum(div), age_category="overall", pop=sum(pop)),
  keyby=eval(grep("age_category", key(tarevents), invert=T, value=T))
]

tarevents <- setnames(setkeyv(rbind(tarevents, tareventsoverall), key(tarevents)), "age_category", "age")

propevents <- tarevents[
  div!=0 | (div==0 & averted==0), {
    vs = averted/div
    vs[is.na(vs)] <- 0
    mn = mean(vs)
    se = 2*sd(vs)/sqrt(.N)
    list(value=mn, CI_low=mn-se, CI_high=mn+se)
  },
  keyby=eval(grep("particle_id", key(tarevents), invert=T, value=T))
][, group:="UF" ][, outcome_denominator := "proportion averted"]

percapevents <- tarevents[, {
    mn = mean(averted/pop*1e5)
    se = 2*sd(averted/pop*1e5)/sqrt(.N)
    list(value=mn, CI_low=mn-se, CI_high=mn+se)
  },
  keyby=eval(grep("particle_id", key(tarevents), invert=T, value=T))
][, group:="UF" ][, outcome_denominator := "per 100,000 pop at risk"]

saveRDS(rbind(propevents,percapevents),paste0(dbpath,"/CMDVI/Phase II analysis/Data/UF-Longini/longini-vacimpact.rds"))

# rm(ls=list())
# cbPalette <- c("Hopkins/UF"="#999999",
#                "Imperial"="#E69F00",
#                "Duke"="#56B4E9",
#                "Notre Dame"="#009E73",
#                "UF"="#F0E442",
#                "Exeter/Oxford"="#0072B2",
#                "Sanofi Pasteur"="#D55E00",
#                "UWA"="#CC79A7",
#                "A"="white",
#                "Longini-ODE"="black"); 
# 
# noGroup_cbPalette <- c("#999999","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","black")
# 
# scen="reference"
# df=subset(rbind(propevents,percapevents), scenario==scen )
# 
# #vaccine impact
# for (my_year in c("cum10","cum30")){
#   #my_outcome_denominator <- "proportion averted"
#   for (my_outcome_denominator in c("proportion averted","per 100,000 pop at risk")){
#     
#     df_tmp=subset(df, (age %in% c("overall", "<9yrs", "9-18yrs", "19+yrs")) & year==my_year &  
#                     outcome_denominator==my_outcome_denominator)
#     dodge=position_dodge(width=.6)
#     
#     p=ggplot(df_tmp, aes(x=age, y=value, ymin=CI_low, ymax=CI_high, fill=group, color=group)) +
#       geom_bar(stat="identity", position=dodge, alpha=0.5, color=NA) +
#       geom_errorbar(stat="identity", position=dodge, alpha=0.8, width=0, size=0.33) +
#       facet_grid(outcome~transmission_setting, scale="free_y") +
#       xlab("age group") +
#       ylab(my_outcome_denominator) +
#       theme_bw() +
#       theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
#       scale_fill_manual(values=cbPalette) + scale_color_manual(values=cbPalette) + ggtitle(my_year)
#     
#     if(my_outcome_denominator=="proportion averted") p <- p+ coord_cartesian( ylim =c(-0.5,0.5))
#     
#     print(p)
#   }
# }
