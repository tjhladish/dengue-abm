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

tar <- "~/Downloads/who-feb-2016/who-feb-2016-aggregated/"
setwd(tar)

poppath <- "~/git/dengue/pop-merida/pop-merida/population-merida.txt" # needs to be merida instead
seropath <- "~/Downloads/who-feb-2016/auto_output/" # path to serostatus processed logs - run process

age_cats <- c("<9yrs","9-18yrs","19+yrs","overall")

pop <- fread(poppath)[,list(pop=.N),keyby=age][,
  age_category := factor(ifelse(
    age < 9, age_cats[1], ifelse(
      age < 19, age_cats[2],
      age_cats[3]
    )), levels = age_cats, ordered = TRUE)
]

vacfiles <- list.files(pattern = "^[0-4][01]{2}1[01]{3}\\.")

bg9files <- list.files(pattern = "^[0-4]0+\\.")

seedkey <- setkey(rbindlist(lapply(vacfiles, function(nm) {
  c(list(seed=as.numeric(gsub(".+\\.","", nm))), fread(nm, nrows = 1, colClasses = c(scenario="character"))[,list(serial, scenario)])
})), seed)[, particle_id := floor(serial/80) ]

cohortdata <- setkey(subset(
  setkey(fread("../cohort.out", col.names = c("seed","year","serostatus","vaccination","severity","count")), seed, serostatus, severity, vaccination)[seedkey],
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

cohortdata <- translate_scenario(cohortdata)

cohortdata[severity!=0, infections := count]
cohortdata[severity>1, symptomatic_cases := count ]
cohortdata[severity>2, hospitalised_cases := count ]
cohortdata[severity==2, c("hospitalised_cases","deaths") := {
  h = rbinom(.N, count, phc)
  d = rbinom(.N, h, cfr)
  list(h,d)
}]
cohortdata[severity==3, deaths := rbinom(.N, hospitalised_cases, cfr)]

overpop <- cohortdata[,list(pop=sum(count)),keyby=list(particle_id, scenario, transmission_setting, year, vaccination, serostatus)]
cohortres <- cohortdata[overpop][severity != 0,
  list(
    pop100k=pop[1]/1e5,
    infections=sum(infections), symptomatic_cases=sum(symptomatic_cases),
    hospitalised_cases=sum(hospitalised_cases), deaths=sum(deaths)
  ), keyby = list(scenario, transmission_setting, particle_id, year, serostatus, vaccination)
]

vacpop <- cohortdata[vaccination==1,list(vacpop=sum(count)),keyby=list(particle_id, scenario, transmission_setting, year)][overpop][,list(vac_percent=vacpop/pop),keyby=list(particle_id, scenario, transmission_setting, year)]

overallcohort <- cohortres[,
  list(
    infections=sum(infections), symptomatic_cases=sum(symptomatic_cases),
    hospitalised_cases=sum(hospitalised_cases), deaths=sum(deaths), pop100k=sum(pop100k)),
  keyby=list(particle_id, transmission_setting, year, scenario)
]

fakeoverall <- setkey(overallcohort[year == 0, list(infections, symptomatic_cases, hospitalised_cases, deaths) ,keyby=key(overallcohort)][,{
  infs  = sample(infections, 30, rep=T)
  symps = rbinom()
  hosps = rbinom(30, symps, 0.5)
  dea = rbinom(30, hosps, cfr)
  list(year=0:29, infections=infs, symptomatic_cases=sample(symptomatic_cases, 30, rep=T),
       hospitalised_cases=sample(hospitalised_cases, 30, rep=T), deaths=sample(deaths, 30, rep=T)) 
}, keyby=list(particle_id, transmission_setting)
], particle_id, transmission_setting, year)

overallcohortres <- overallcohort[fakeoverall][,
  list(
    infections_averted = (i.infections - infections)/pop100k,
    symptomatic_cases_averted = (i.symptomatic_cases - symptomatic_cases)/pop100k,
    hospitalised_cases_averted = (i.hospitalised_cases - hospitalised_cases)/pop100k,
    deaths_averted = (i.deaths - deaths)/pop100k
  ),
  keyby=key(overallcohort)
][,
  list(year,
    infections_averted = cumsum(infections_averted),
    symptomatic_cases_averted = cumsum(symptomatic_cases_averted),
    hospitalised_cases_averted = cumsum(hospitalised_cases_averted),
    deaths_averted = cumsum(deaths_averted)
  ),
  keyby=list(particle_id, transmission_setting, scenario)
]


melt(, id.vars = key(overallcohort), variable.name = "outcome")[,{
  
    list(age="cohort_overall", outcome_denominator="cumulative - per 100,000 pop at risk")
  },
  keyby=list(scenario, transmission_setting, year, outcome)
]



cohortslice <- function(dt, agelab, ...) subset(dt[...][, age:=agelab ], select=-c(serostatus, vaccination))

negun <- cohortslice(cohortres, "cohort_vaccneg_seroneg", serostatus==0 & vaccination==0)
posun <- cohortslice(cohortres, "cohort_vaccneg_seropos", serostatus==1 & vaccination==0)
negvac <- cohortslice(cohortres, "cohort_vaccpos_seroneg",serostatus==0 & vaccination==1)
posvac <- cohortslice(cohortres, "cohort_vaccpos_seropos",serostatus==1 & vaccination==1)
## TODO replace this with real data
# chort <- cohortres[
#   year == 0, count, keyby=list(scenario, transmission_setting, particle_id, serostatus, severity, vaccination)
# ][,
#   list(count=sum(count)), keyby=list(scenario, transmission_setting, particle_id, serostatus, severity)
# ][,
#   list(count, year=0:29), keyby=key(chort)
# ]


cohortres[,list(year, cumsum(count)/cumsum(pop)),keyby=list(scenario, particle_id, serostatus, vaccination, severity)]


econdata <- subset(merge(
  dcast.data.table(melt(
    rbindlist(lapply(c(vacfiles, bg9files), fread, colClasses = c(scenario="character"))),
      id.vars = c("serial","scenario","age","outcome"), variable.name = "year", value.name = "count"
    )[outcome != 0][,
      year := as.integer(gsub("y", "", year))
    ][,
      outcome := factor(c("mild","severe")[outcome], levels=c("mild","severe"), ordered = T)
    ], serial + scenario + age + year ~ outcome, value.var = "count"
  ), pop, by = "age"), select=-age_category)[,
  particle_id := floor(serial/80)
][,
  symptomatic := round(mild*(1-phc))
][,
  amb := symptomatic
][,
  hosp := round(((mild-amb)+severe)*(1-cfr))
][,
  death := (mild + severe) - (amb+hosp)
][,
  transmission_setting := seq(10,90,20)[as.integer(substr(scenario, 1, 1))+1]
][,
  vacc_coverage := ifelse(substr(scenario,4,4) == 0, 0, ifelse(substr(scenario,7,7) == 0, .5, .8))
]



econdata[vacc_coverage == 0.0, scen := "noVaccine"]
econdata[vacc_coverage == 0.5, scen := "coverageAT50%"]
econdata[substr(scenario,3,3) == "1", scen := "altVaccine" ]
econdata[substr(scenario,2,2) == "1", scen := "catchUp"]
econdata[substr(scenario,6,6) == "1", scen := paste0(scen,"To30")]
econdata[grepl("^[0-4]0{2}10{2}1",scenario), scen := "reference"]
econdata[substr(scenario,5,5) == "1", scen := "routineAT16"]

econdata[, vac := 0]
econdata[age==9 & scen != "routineAT16", vac := round(vacc_coverage*pop)]
econdata[age %in% (10:17) & scen == "catchUp" & year == 0, vac := round(vacc_coverage*pop)]
econdata[age %in% (10:30) & scen == "catchUpTo30" & year == 0, vac := round(vacc_coverage*pop)]
econdata[age==16 & scen == "routineAT16", vac := round(vacc_coverage*pop)]

econfinal <- subset(econdata, select=c(particle_id, scen, transmission_setting, year, age, vac, amb, hosp, death))
setkey(econfinal, particle_id, scen, transmission_setting, year, age)

saveRDS(econfinal, "~/Downloads/econref.rds")

bg16files <- list.files(pattern = "^[01]0+\\.")
vac9files <- list.files(pattern = "^[0-4][01]{2}10[01]{2}\\.")
vac16files <- list.files(pattern = "^[01][01]{2}11[01]{2}\\.")

pop_ref <- pop[,list(pop=sum(pop)/1e5),by=age_category]
pop_ref <- rbind(pop_ref, pop_ref[,list(pop=sum(pop), age_category="overall")])

# cohort1 <- pop[age >= 9][1:30, list(pop, year=0:29, age="cohort1")]
# cohort2 <- pop[age >= 9][1:29, list(pop, year=1:29, age="cohort2")]
# cohort3 <- pop[age >= 9][1:28, list(pop, year=2:29, age="cohort3")]
# cohort4 <- pop[age >= 9][1:27, list(pop, year=3:29, age="cohort4")]
# cohort5 <- pop[age >= 9][1:26, list(pop, year=4:29, age="cohort5")]

# cohorts <- cohort1 # rbind(cohort1, cohort2, cohort3, cohort4, cohort5)
# cohorts[, year:=factor(year)][,pop:=pop/1e5]

## TODO get sample size, divide CIs by sqrt sample size == should be basically /10

munger <- function(fname, vacage=9, xmis_setting=seq(10,90,by=20)) {
  long <- subset(melt(fread(fname, colClasses = c(scenario="character")),
   id.vars = c("serial","scenario","age","outcome"),
   variable.name = "year", value.name = "count")[,
   year := as.integer(gsub("y","",year))
  ][,
    event := factor(
      c("asymptomatic","symptomatic","severe")[outcome+1],
      levels = c("asymptomatic","symptomatic","severe"), ordered = T
    )
  ][,
    particle_id := floor(serial/80)
  ][,
    transmission_setting := xmis_setting[as.integer(substr(fname,1,1))+1]
  ], select = -serial)[,
    scenario := substr(scenario,2,7)
  ]
  
  imp_breakdown <- long[,
    age_category := factor(ifelse(
      age < 9, age_cats[1], ifelse(
        age < 19, age_cats[2],
        age_cats[3]
      )), levels = age_cats, ordered = TRUE)
    ][,list(count=sum(count)),by=list(particle_id, transmission_setting, scenario, year, event, age=age_category)]
  
  imp_overall <- imp_breakdown[,list(age = factor(age_cats[4],levels=age_cats,ordered = T), count=sum(count)),by=list(particle_id, transmission_setting, scenario, year, event)]
  
  impact <- rbind(imp_breakdown, imp_overall)
  # ref$impact - impact == cases averted
  # del / ref$impact == proportion
#   long[,cohort:="none"]
#   long[age-year == vacage, cohort := "cohort1"]
#   long[age-year == vacage-1 & year > 0, cohort := "cohort2"]
#   long[age-year == vacage-2 & year > 1, cohort := "cohort3"]
#   long[age-year == vacage-3 & year > 2, cohort := "cohort4"]
#   long[age-year == vacage-4 & year > 3, cohort := "cohort5"]
#   cohorts <- long[cohort != "none",count,by=list(particle_id, transmission_setting, scenario,year,event,age=cohort)]
  
  rbindlist(list(
    overall_allyears = impact[age=="overall", ],
    cum1030=setcolorder(melt(impact[,
        list(cum10=sum(count[year<10]), cum30=sum(count)),
        by = list(particle_id, transmission_setting, scenario, event, age)
      ], measure.vars = c("cum10","cum30"), variable.name = "year", value.name = "count"),
      c("particle_id", "transmission_setting", "scenario","year","event","age","count")) #,
    #cohorts = cohorts
  ))
}

backgroundRead <- function(files, va=9) setkey(
  rbindlist(mclapply(
    files, munger, vacage=va, mc.cores = detectCores()-1
  )), particle_id, transmission_setting, year, event, age, scenario)

bg9 <- subset(backgroundRead(bg9files), select=-scenario)
bg16 <- subset(backgroundRead(bg16files, va=16), select=-scenario)

vac9 <- backgroundRead(vac9files)
vac16 <- backgroundRead(vac16files, va=16)

items <- rbind(vac9[bg9], vac16[bg16])[
  ## convert asymptomatics to infections
  ## convert symptomatics -> infections vs hospitalizations
  ## convert severe -> hospitalizations
][,
  list(averted=i.count-count, vs=i.count),
  by=list(particle_id, transmission_setting, scenario, year, event, age)
]
class(items$vs) <- "double"
items[vs == 0 & !(age %in% c("cohort1","cohort2","cohort3","cohort4","cohort5")), vs := 1e-1]

hosp <- items[event != "asymptomatic",
  list(
    outcome="hospitalised_cases_averted",
    averted=sum(averted*ifelse(event == "symptomatic", phc, 1)),
    vs=sum(vs*ifelse(event == "symptomatic", phc, 1))
  ),
  by=list(particle_id, transmission_setting, scenario, year, age)
]
death <- items[year %in% c("cum10","cum30") & event != "asymptomatic",
  list(
    outcome="deaths_averted",
    averted=sum(averted)*cfr,
    vs=sum(vs)*cfr
  ),
  by=list(particle_id, transmission_setting, scenario, year, age)
]

newitems <- rbind(
  items[, list(outcome="infections_averted", averted=sum(averted), vs=sum(vs)), by=list(particle_id, transmission_setting, scenario, year, age)],
  items[event != "asymptomatic", list(outcome="symptomatic_cases_averted", averted=sum(averted), vs=sum(vs)), list(particle_id, transmission_setting, scenario, year, age)],
  hosp,
  death
)

newitems[year %in% c("cum10","cum30"), vs := averted/vs ]
newitems[!(year %in% c("cum10","cum30")), vs := NA ]

moltenitems <- melt(newitems, measure.vars = c("averted", "vs"))[!is.na(value)]
nonsero <- rbind(
  subset(merge(moltenitems[variable=="averted"], pop_ref[,pop,keyby=list(age=age_category)], by="age")[, value := value/pop ], select=-pop),
  #subset(merge(moltenitems[variable=="averted"], cohorts, by=c("age","year"))[, value := value/pop ], select=-pop),
  setcolorder(moltenitems[variable != "averted"], c("age", "particle_id", "transmission_setting", "scenario", "year", "outcome", "variable", "value"))
)
# reduce population here

reduceditems <- nonsero[,{
  v = mean(value)
  sdv = sd(value)
  list(value=v,CI_low=v-sdv,CI_high=v+sdv)
},by=list(transmission_setting, scenario, year, outcome, age, variable)]

reduceditems[, outcome_denominator := ifelse(variable == "averted", "per 100,000 pop at risk", "proportion averted") ]

translate_scenario(reduceditems)

reduceditems[, group:= "longini" ]
reduceditems<-subset(reduceditems, select=-c(variable, scen))
require(bit64)
seroscn <- setkey(subset(fread(paste0(seropath,"seroscenarios.csv")), select=c(V1, V9, V10, V11, V12, V13, V14, V15))[,
  unique(V1), keyby=list(V9,V10,V11,V12,V13,V14,V15)
][, transmission_setting := seq(10,90,by=20)[V9+1] ][,
  scenario := paste0(V10, V11, V12, V13, V14, V15)
][,
  list(seed=unique(V1)), by=list(transmission_setting, scenario)
], seed)

seroscn[substr(scenario,6,6) == "0", scen := "coverageAT50%"]
seroscn[substr(scenario,2,2) == "1", scen := "altVaccine" ]
seroscn[substr(scenario,1,1) == "1", scen := "catchUp"]
seroscn[substr(scenario,5,5) == "1", scen := paste0(scen,"To30")]
seroscn[grepl("^0{2}10{2}1",scenario), scen := "reference"]
seroscn[substr(scenario,4,4) == "1", scen := "routineAT16"]
seroscn[substr(scenario,3,3) == "0", scen := "noVaccine"]
seroscn[, scenario:=scen ]

serosrv <- setkey(fread(paste0(seropath,"seroprevalence.csv"))[,
  list(seed=unique(seed)), keyby=list(year, vax_age, seroprevalence)
][, year := year-50 ], seed, year)

seroresults <- setcolorder(
  serosrv[seroscn][, 
    list(value = mean(seroprevalence), CI_low=NA, CI_high=NA, outcome="seropositive", outcome_denominator="proportion", group="longini"),
    by=list(age=paste0(vax_age,"yrs"), transmission_setting, scenario, year)
  ],
  names(reduceditems)
)

## TODO using baseline case data, determine proportion of symptomatic + hospitalized cases by age, averaged over...saaaay...10 years?
##  then turn that into cumulative cases by those ages
##  outcome = symptomatic | hospitalised cases, outcome_denominator == cumulative proportion at baseline?

## TODO do cumulative X_averted for everything, rbind it, have outcome_denom = "cumulative - per 100,000 pop at risk"

saveRDS(rbind(reduceditems, seroresults), "~/Dropbox/CMDVI/Phase II analysis/Data/UF-Longini/longini2.RData")

# mort_rate <- c(symptomatic=0, severe=0.1)
# hosp_rate <- c(symptomatic=0.111)
# hosp|symp = hosp|symp,pri + symp,sec = 
# hosp|severe = hosp|sev,pri + sev,sec
# hosp|pri = hosp | pri,symp + pri,severe
# hosp|sec = 1.88*hosp|pri = hosp | sec,symp + sec,severe

## some checks
# bg <- melt(bg9[year %in% 0:29 & age == "overall",{
#   res <- as.list(quantile(count, probs = c(.1,.5,.9), names=F))
#   names(res) <- c("low","med","hi")
#   res
# }, by=list(transmission_setting, event, year)], measure.vars = c("low","med","hi"))[variable=="med"][, scn := "background" ]
# fg80 <- melt(vac9[year %in% 0:29 & age == "overall" & scenario == "001001",{
#   res <- as.list(quantile(count, probs = c(.1,.5,.9), names=F))
#   names(res) <- c("low","med","hi")
#   res
# }, by=list(transmission_setting, event, year)], measure.vars = c("low","med","hi"))[variable=="med"][, scn := "foreground_80" ]
# fg80_catchup <- melt(vac9[year %in% 0:29 & age == "overall" & scenario == "101001",{
#   res <- as.list(quantile(count, probs = c(.1,.5,.9), names=F))
#   names(res) <- c("low","med","hi")
#   res
# }, by=list(transmission_setting, event, year)], measure.vars = c("low","med","hi"))[variable=="med"][, scn := "foreground_80_catchup" ]
# 
# fg80_catchup30 <- melt(vac9[year %in% 0:29 & age == "overall" & scenario == "101011",{
#   res <- as.list(quantile(count, probs = c(.1,.5,.9), names=F))
#   names(res) <- c("low","med","hi")
#   res
# }, by=list(transmission_setting, event, year)], measure.vars = c("low","med","hi"))[variable=="med"][, scn := "foreground_80_catchup30" ]
# 
# fg50 <- melt(vac9[year %in% 0:29 & age == "overall" & scenario == "001000",{
#   res <- as.list(quantile(count, probs = c(.1,.5,.9), names=F))
#   names(res) <- c("low","med","hi")
#   res
# }, by=list(transmission_setting, event, year)], measure.vars = c("low","med","hi"))[variable=="med"][, scn := "foreground_50" ]
# 
# ggplot(rbind(bg,fg80, fg80_catchup, fg80_catchup30, fg50)) +
#   theme_bw() +
#   aes(x=as.numeric(year), y=value, color=scn) +
#   facet_grid(event ~ transmission_setting, scales = "free_y") + geom_line() +
#   stat_smooth(method="lm", se = F)
# 
# ggplot(melt(vac9[year %in% 0:29 & age == "overall" & scenario == "001001",{
#   res <- as.list(quantile(count, probs = c(.1,.5,.9), names=F))
#   names(res) <- c("low","med","hi")
#   res
# }, by=list(transmission_setting, event, year)], measure.vars = c("low","med","hi"))) +
#   theme_bw() +
#   aes(x=as.numeric(year), y=value) +
#   facet_grid(event ~ transmission_setting, scales = "free_y") + geom_line()