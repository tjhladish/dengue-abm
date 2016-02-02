tar <- "~/Downloads/who-feb-2016/who-feb-2016-aggregated/"
setwd(tar)

rm(list=ls())
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

poppath <- "/Volumes/Data/workspaces/dengue/pop-merida/pop-merida/population-merida.txt" # needs to be merida instead
seropath <- "~/Downloads/who-feb-2016/auto_output/" # path to serostatus processed logs - run process

age_cats <- c("<9yrs","9-18yrs","19+yrs","overall")

pop <- fread(poppath)[,list(pop=.N),keyby=age][,
  age_category := factor(ifelse(
    age < 9, age_cats[1], ifelse(
      age < 19, age_cats[2],
      age_cats[3]
    )), levels = age_cats, ordered = TRUE)
]

srcfiles <- list.files(pattern = "^[0-4].+\\..+")

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

saveRDS(cached_res, "~/Downloads/who-rawpre.rds")

cumulative_events <- setkeyv(cached_res[,
  list(year=paste0("cum", year+1), value=cumsum(value), i.value=cumsum(i.value)),
  keyby=eval(grep("year", key(cached_res), invert = T, value=T))
][, averted := i.value - value ][, outcome := paste0(outcome, "_averted")], key(cached_res))

tarevents <- cumulative_events[year %in% c("cum10","cum30")][,
  list(div=sum(i.value), averted=sum(averted), pop=sum(pop)),
  keyby=eval(c(grep("(age|pop)",key(cumulative_events),invert=T,value=T),"age_category"))
]

propevents <- tarevents[
  div!=0 | (div==0 & averted==0), {
    vs = averted/div
    vs[is.na(vs)] <- 0
    mn = mean(vs)
    se = 2*sd(vs)/sqrt(.N)
    list(value=mn, CI_low=mn-se, CI_high=mn+se)
  },
  keyby=eval(grep("particle_id", key(tarevents), invert=T, value=T))
]

df=subset(df_all, scenario==scen )

#vaccine impact
for (my_year in c("cum10","cum30")){
  for (my_outcome_denominator in c("proportion averted","per 100,000 pop at risk")){
    
    df_tmp=subset(df, (age %in% c("overall", "<9yrs", "9-18yrs", "19+yrs")) & year==my_year &  
                    outcome_denominator==my_outcome_denominator)
    dodge=position_dodge(width=.6)
    
    p=ggplot(df_tmp, aes(x=age, y=value, ymin=CI_low, ymax=CI_high, fill=group, color=group)) +
      geom_bar(stat="identity", position=dodge, alpha=0.5, color=NA) +
      geom_errorbar(stat="identity", position=dodge, alpha=0.8, width=0, size=0.33) +
      facet_grid(outcome~transmission_setting, scale="free_y") +
      xlab("age group") +
      ylab(my_outcome_denominator) +
      theme_bw() +
      theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
      scale_fill_manual(values=cbPalette) + scale_color_manual(values=cbPalette)
    
    if(my_outcome_denominator=="proportion averted") p <- p+ coord_cartesian( ylim =c(-0.5,0.5))
    
    p
  }
}


# backgroundRead <- function(files, va=9) setkey(
#   rbindlist(mclapply(
#     files, munger, vacage=va, mc.cores = detectCores()-1
#   )), particle_id, transmission_setting, year, event, age, scenario)
# 
# bg9 <- subset(backgroundRead(bg9files), select=-scenario)
# bg16 <- subset(backgroundRead(bg16files, va=16), select=-scenario)
# 
# vac9 <- backgroundRead(vac9files)
# vac16 <- backgroundRead(vac16files, va=16)
# 
# items <- rbind(vac9[bg9], vac16[bg16])[
#   ## convert asymptomatics to infections
#   ## convert symptomatics -> infections vs hospitalizations
#   ## convert severe -> hospitalizations
# ][,
#   list(averted=i.count-count, vs=i.count),
#   by=list(particle_id, transmission_setting, scenario, year, event, age)
# ]
# class(items$vs) <- "double"
# items[vs == 0 & !(age %in% c("cohort1","cohort2","cohort3","cohort4","cohort5")), vs := 1e-1]
# 
# hosp <- items[event != "asymptomatic",
#   list(
#     outcome="hospitalised_cases_averted",
#     averted=sum(averted*ifelse(event == "symptomatic", phc, 1)),
#     vs=sum(vs*ifelse(event == "symptomatic", phc, 1))
#   ),
#   by=list(particle_id, transmission_setting, scenario, year, age)
# ]
# death <- items[year %in% c("cum10","cum30") & event != "asymptomatic",
#   list(
#     outcome="deaths_averted",
#     averted=sum(averted)*cfr,
#     vs=sum(vs)*cfr
#   ),
#   by=list(particle_id, transmission_setting, scenario, year, age)
# ]
# 
# newitems <- rbind(
#   items[, list(outcome="infections_averted", averted=sum(averted), vs=sum(vs)), by=list(particle_id, transmission_setting, scenario, year, age)],
#   items[event != "asymptomatic", list(outcome="symptomatic_cases_averted", averted=sum(averted), vs=sum(vs)), list(particle_id, transmission_setting, scenario, year, age)],
#   hosp,
#   death
# )
# 
# newitems[year %in% c("cum10","cum30"), vs := averted/vs ]
# newitems[!(year %in% c("cum10","cum30")), vs := NA ]
# 
# moltenitems <- melt(newitems, measure.vars = c("averted", "vs"))[!is.na(value)]
# nonsero <- rbind(
#   subset(merge(moltenitems[variable=="averted"], pop_ref[,pop,keyby=list(age=age_category)], by="age")[, value := value/pop ], select=-pop),
#   #subset(merge(moltenitems[variable=="averted"], cohorts, by=c("age","year"))[, value := value/pop ], select=-pop),
#   setcolorder(moltenitems[variable != "averted"], c("age", "particle_id", "transmission_setting", "scenario", "year", "outcome", "variable", "value"))
# )
# # reduce population here
# 
# reduceditems <- nonsero[,{
#   v = mean(value)
#   sdv = sd(value)
#   list(value=v,CI_low=v-sdv,CI_high=v+sdv)
# },by=list(transmission_setting, scenario, year, outcome, age, variable)]
# 
# reduceditems[, outcome_denominator := ifelse(variable == "averted", "per 100,000 pop at risk", "proportion averted") ]
# 
# translate_scenario(reduceditems)
# 
# reduceditems[, group:= "longini" ]
# reduceditems<-subset(reduceditems, select=-c(variable, scen))
# require(bit64)
# seroscn <- setkey(subset(fread(paste0(seropath,"seroscenarios.csv")), select=c(V1, V9, V10, V11, V12, V13, V14, V15))[,
#   unique(V1), keyby=list(V9,V10,V11,V12,V13,V14,V15)
# ][, transmission_setting := seq(10,90,by=20)[V9+1] ][,
#   scenario := paste0(V10, V11, V12, V13, V14, V15)
# ][,
#   list(seed=unique(V1)), by=list(transmission_setting, scenario)
# ], seed)
# 
# seroscn[substr(scenario,6,6) == "0", scen := "coverageAT50%"]
# seroscn[substr(scenario,2,2) == "1", scen := "altVaccine" ]
# seroscn[substr(scenario,1,1) == "1", scen := "catchUp"]
# seroscn[substr(scenario,5,5) == "1", scen := paste0(scen,"To30")]
# seroscn[grepl("^0{2}10{2}1",scenario), scen := "reference"]
# seroscn[substr(scenario,4,4) == "1", scen := "routineAT16"]
# seroscn[substr(scenario,3,3) == "0", scen := "noVaccine"]
# seroscn[, scenario:=scen ]
# 
# serosrv <- setkey(fread(paste0(seropath,"seroprevalence.csv"))[,
#   list(seed=unique(seed)), keyby=list(year, vax_age, seroprevalence)
# ][, year := year-50 ], seed, year)
# 
# seroresults <- setcolorder(
#   serosrv[seroscn][, 
#     list(value = mean(seroprevalence), CI_low=NA, CI_high=NA, outcome="seropositive", outcome_denominator="proportion", group="longini"),
#     by=list(age=paste0(vax_age,"yrs"), transmission_setting, scenario, year)
#   ],
#   names(reduceditems)
# )

## TODO using baseline case data, determine proportion of symptomatic + hospitalized cases by age, averaged over...saaaay...10 years?
##  then turn that into cumulative cases by those ages
##  outcome = symptomatic | hospitalised cases, outcome_denominator == cumulative proportion at baseline?

## TODO do cumulative X_averted for everything, rbind it, have outcome_denom = "cumulative - per 100,000 pop at risk"

#saveRDS(rbind(reduceditems, seroresults), "~/Dropbox/CMDVI/Phase II analysis/Data/UF-Longini/longini2.RData")


