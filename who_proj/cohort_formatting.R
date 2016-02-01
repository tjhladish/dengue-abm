## cohort review

wd <- "~/Dropbox/who-feb-2016/who-feb-2016-aggregated"
setwd(wd)

rm(list=ls())

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

cfr <- .00078

require(reshape2)
require(data.table)
require(bit64)

vacfiles <- list.files(pattern = "^[0-4][01]{2}1[01]{3}\\.")
bgfiles <- list.files(pattern = "^[0-4]0+\\.")

seedkeyFG <- setkey(rbindlist(lapply(vacfiles, function(nm) {
  c(list(seed=as.numeric(gsub(".+\\.","", nm))), fread(nm, nrows = 1, colClasses = c(scenario="character"))[,list(serial, scenario)])
})), seed)[, particle_id := floor(serial/80) ]

seedkeyBG <- setkey(rbindlist(lapply(bgfiles, function(nm) {
  c(list(seed=as.numeric(gsub(".+\\.","", nm))), fread(nm, nrows = 1, colClasses = c(scenario="character"))[,list(serial, scenario)])
})), seed)[, particle_id := floor(serial/80) ]

readInCohort <- function(src, sk, cnames = c("seed","year","serostatus","vaccination","severity","count")) {
  res <- setkey(subset(
    setkey(fread(src, col.names = cnames), seed, serostatus, severity, vaccination)[sk],
    select = -c(seed, serial)
  )[,
    `:=`(
      transmission_setting = seq(10,90,20)[as.integer(substr(scenario,1,1))+1],
      scenario = substr(scenario,2,7)
    )                                                                 
    ], particle_id, scenario, transmission_setting, year, vaccination, serostatus, severity)[,
   `:=`(
      infections = 0,
      symptomatic_cases = 0,
      hospitalised_cases = 0,
      deaths = 0
    )                                                                                       
  ]
  class(res$count) <- "double"
  res
}

cohortdata <- readInCohort("../cohort/cohort.vaccinations", seedkeyFG)
refcohortdata <- readInCohort(
  "../cohort/cohort.no-vaccinations", seedkeyBG,
  c("seed","year","serostatus","vaccination","severity","count","routineage")
)[vaccination == 0]

process_severity <- function(tar) {
  tar[severity!=0, infections := count] # 0 == non infection
  tar[severity>1, symptomatic_cases := count ] # 1 == asymptomatic infection
  tar[severity==2, hospitalised_cases := count*phc ] # 2 == mild disease, possible hosp
  tar[severity==3, hospitalised_cases := count ] # 3 == severe disease, assume hosp
  tar[severity>=2, deaths := hospitalised_cases*cfr] # death is % of hosp
}

cohortify <- function(src) {
  pop <- src[,
    list(pop=sum(count)),
    keyby=list(scenario, transmission_setting, particle_id, year, serostatus, vaccination)
  ]
  melt(src,
    measure.vars = c("infections", "symptomatic_cases", "hospitalised_cases", "deaths"),
    variable.name = "outcome"
  )[,
    list(value=sum(value)),
    keyby=list(scenario, transmission_setting, particle_id, year, serostatus, vaccination, outcome)
  ][pop]
}

fincohortdata <- cohortify(translate_scenario(process_severity(cohortdata)))
refcohortdata <- translate_scenario(process_severity(refcohortdata))

## duplicate refcohort data, with vaccine == 1
repcohdata <- copy(refcohortdata)[, vaccination := 1 ]

scen16 <- "routineAT16"
nonscen16 <- fincohortdata[scenario != scen16,unique(scenario)]

tot <- setkeyv(rbind(refcohortdata, repcohdata), key(refcohortdata))
regref <- tot[routineage==9]

regrefexp <- regref[,
  list(scenario=nonscen16, infections, symptomatic_cases, hospitalised_cases, deaths, count),
  keyby=eval(grep("scenario", key(regref), value=T, invert=T))
]
ref16 <- subset(tot[routineage==16][, scenario := scen16 ], select=-routineage)
setcolorder(regrefexp, names(ref16))
refcoho <- cohortify(setkeyv(rbind(regrefexp, ref16), key(refcohortdata)))

averted <- fincohortdata[refcoho][,
  list(value = (ifelse(i.pop==0,0,i.value/i.pop) - ifelse(pop==0,0,value/pop))*1e5),
  keyby=key(fincohortdata)
]

cumaverted <- averted[,
  list(value=cumsum(value), year), keyby=eval(grep("year", key(averted), value=T, invert=T))
]

totstatcumaverted <- cumaverted[,
  list(value=sum(value)), keyby=list(scenario, transmission_setting, particle_id, year, outcome)                            
][,{
    mn = mean(value)
    se = 2*sd(value)/sqrt(.N)
    list(value=mn, CI_low = mn-se, CI_high = mn + se)
  }, keyby=list(scenario, transmission_setting, year, outcome)
][, age:= "cohort_complete" ]

statcumaverted <- cumaverted[,{
    mn = mean(value)
    se = 2*sd(value)/sqrt(.N)
    list(value=mn, CI_low = mn-se, CI_high = mn + se)
  },
  keyby=list(scenario, transmission_setting, year, serostatus, vaccination, outcome)
]

statcumaverted[,
  serostatus:= ifelse(serostatus==1,"pos","neg")
][,
  vaccination:= ifelse(vaccination==1,"pos","neg")
][,
  age := paste0("cohort_vacc",vaccination,"_sero",serostatus)
]

output <- rbind(
  statcumaverted[(year == 9) | (year == 29)][,
    list(value, CI_low, CI_high),
    keyby=list(scenario, transmission_setting, year=paste0("cum",year+1), outcome, age)
  ],
  totstatcumaverted
)[, outcome_denominator:="cumulative - per 100,000 pop at risk"][, group:="longini"][, outcome := paste0(outcome,"_averted")]

saveRDS(output, "~/Dropbox/CMDVI/Phase II analysis/Data/UF-Longini/longini-cohort.rds")

# my_age= "cohort_complete"; cbPalette <- c("red","green","blue")
# df_tmp=subset(output, age == my_age & (outcome_denominator=="cumulative - per 100,000 pop at risk") & scenario == "reference" & year %in% 1:30 )
# df_tmp$year=as.numeric(as.character(df_tmp$year))
# 
# ggplot(df_tmp, aes(x= year, y=value, ymin=CI_low, ymax=CI_high, fill=group, color=group, group=group)) +
#   geom_line() +
#   geom_ribbon(alpha=0.2, color=NA) +
#   facet_grid(outcome~transmission_setting, scale="free_y") +
#   xlab("time after introduction of CYD (years)") +
#   ylab("cumulative - per 100,000 pop at risk") +
#   theme_bw() +
#   scale_fill_manual(values=cbPalette) + scale_color_manual(values=cbPalette) 
# 
# #vaccine effects in vaccinated cohort split up
# df_tmp=subset(output, (age %in% c("cohort_complete","cohort_vaccpos_seropos","cohort_vaccpos_seroneg","cohort_vaccneg_seropos","cohort_vaccneg_seroneg")) & 
#                 (outcome_denominator=="cumulative - per 100,000 pop at risk") & (year %in% c("cum10","cum30") & scenario == "reference") )
# df_tmp$age=gsub("cohort_complete","complete\ncohort",df_tmp$age)
# df_tmp$age=gsub("cohort_vaccpos_seropos","vaccinated\nseropositive",df_tmp$age)
# df_tmp$age=gsub("cohort_vaccpos_seroneg","vaccinated\nseronegative",df_tmp$age)
# df_tmp$age=gsub("cohort_vaccneg_seropos","unvaccinated\nseropositive",df_tmp$age)
# df_tmp$age=gsub("cohort_vaccneg_seroneg","unvaccinated\nseronegative",df_tmp$age)
# ggplot(df_tmp, aes(x=age, y=value, ymin=CI_low, ymax=CI_high, fill=group, color=group, shape=year, linetype=year)) +
#   geom_pointrange(position=position_dodge(.5)) +
#   facet_grid(outcome~transmission_setting, scale="free_y") +
#   xlab("first vaccine eligible cohort") +
#   ylab("per 100,000 pop at risk") +
#   theme_bw() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
#   scale_fill_manual(values=cbPalette) + scale_color_manual(values=cbPalette) 