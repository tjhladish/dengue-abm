rm(list=ls())
require(reshape2)
require(data.table)
require(parallel)
# require(ggplot2)

tar <- "~/Downloads/who-jan-2016-aggregated/"
setwd(tar)

poppath <- "~/git/dengue/pop-merida/pop-merida/population-merida.txt" # needs to be merida instead
seropath <- "~/Downloads/auto_output/" # path to serostatus processed logs - run process

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

econdata <- subset(merge(
  dcast.data.table(melt(
    rbindlist(lapply(c(vacfiles, bg9files), fread, colClasses = c(scenario="character"))),
      id.vars = c("serial","scenario","age","outcome"), variable.name = "year", value.name = "count"
    )[outcome != 0][,
      year := as.integer(gsub("y", "", year))
    ][,
      outcome := factor(c("mild","severe")[outcome], levels=c("mild","severe"), ordered = T)
    ], serial + scenario + age + year ~ outcome, value.var = "count"
  ), pop, by = "age"), select=-age_category)

bg16files <- list.files(pattern = "^[01]0+\\.")
vac9files <- list.files(pattern = "^[0-4][01]{2}10[01]{2}\\.")
vac16files <- list.files(pattern = "^[01][01]{2}11[01]{2}\\.")

pop_ref <- pop[,list(pop=sum(pop)/1e5),by=age_category]
pop_ref <- rbind(pop_ref, pop_ref[,list(pop=sum(pop), age_category="overall")])

cohort1 <- pop[age >= 9][1:30, list(pop, year=0:29, age="cohort1")]
cohort2 <- pop[age >= 9][1:29, list(pop, year=1:29, age="cohort2")]
cohort3 <- pop[age >= 9][1:28, list(pop, year=2:29, age="cohort3")]
cohort4 <- pop[age >= 9][1:27, list(pop, year=3:29, age="cohort4")]
cohort5 <- pop[age >= 9][1:26, list(pop, year=4:29, age="cohort5")]

cohorts <- rbind(cohort1, cohort2, cohort3, cohort4, cohort5)
cohorts[, year:=factor(year)][,pop:=pop/1e5]

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
  long[,cohort:="none"]
  long[age-year == vacage, cohort := "cohort1"]
  long[age-year == vacage-1 & year > 0, cohort := "cohort2"]
  long[age-year == vacage-2 & year > 1, cohort := "cohort3"]
  long[age-year == vacage-3 & year > 2, cohort := "cohort4"]
  long[age-year == vacage-4 & year > 3, cohort := "cohort5"]
  cohorts <- long[cohort != "none",count,by=list(particle_id, transmission_setting, scenario,year,event,age=cohort)]
  
  rbindlist(list(
    overall_allyears = impact[age=="overall", ],
    cum1030=setcolorder(melt(impact[,
      list(cum10=sum(count[year<10]), cum30=sum(count)),
      by = list(particle_id, transmission_setting, scenario, event, age)
    ], measure.vars = c("cum10","cum30"), variable.name = "year", value.name = "count"), c("particle_id", "transmission_setting", "scenario","year","event","age","count")),
    cohorts = cohorts
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

scn_render <- function(scnstr) {
  vacmech <- ifelse(substr(scnstr,2,2) == "0", "baseline", "altvac")
  catchup <- ifelse(substr(scnstr,1,1) == "1", paste0("catchup=>", ifelse(substr(scnstr,5,5) == "0", 17, 30)), "")
  routine <- ifelse(substr(scnstr,4,4) == "0", "", "routine=16yo")
  vaxrate <- ifelse(substr(scnstr,6,6) == "0", "coverage=50", "")
  gsub(",$","",gsub(",+",",",paste(vacmech, catchup, routine, vaxrate, sep=",")))
}

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

# rr_severity_sec_vs_pri = 20
# 
# ph1 = 0.111
# ph2 = ph1*1.88
# ps2 = 20*ps1
# ph1 = phs*ps1 + phc*pc1
# ph2 = phs*ps2 + phc*pc2
# 
# ph1 = 0.111
# ph2 = 0.111*1.88
# ps2 = 20*ps1
# 
# 0.111 = phs*ps1 + phc*(1-ps1)
# 0.111*1.88 = phs*20*ps1 + phc*(1-20*ps1)

# assume phs ~ 1
# (0.111 - ps1)/(1-ps1) = phc
# 0.111*1.88 = 20*ps1 + phc*(1-20*ps1)
# 0.111*(1-1.88) = -19*ps1 - 0.111*1.88*ps1 + 20*0.111*ps1
ps1 = 0.005749711
# ps2 = 20*ps1 = 0.1149942
# pc1 = 1 - ps1
# pc2 = 1 - ps2
phc = (0.111 - ps1)/(1-ps1)

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
    averted=sum(averted)*.00078,
    vs=sum(vs)*.00078
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
  subset(merge(moltenitems[variable=="averted"], cohorts, by=c("age","year"))[, value := value/pop ], select=-pop),
  setcolorder(moltenitems[variable != "averted"], c("age", "particle_id", "transmission_setting", "scenario", "year", "outcome", "variable", "value"))
)
# reduce population here

reduceditems <- nonsero[,{
  v = mean(value)
  sdv = sd(value)
  list(value=v,CI_low=v-sdv,CI_high=v+sdv)
},by=list(transmission_setting, scenario, year, outcome, age, variable)]

reduceditems[, outcome_denominator := ifelse(variable == "averted", "per 100,000 pop at risk", "proportion averted") ]
reduceditems[, scenario:=scn_render(scenario)]

reduceditems[, group:= "longini" ]
reduceditems<-subset(reduceditems, select=-variable)
require(bit64)
seroscn <- setkey(subset(fread(paste0(seropath,"seroscenarios.csv")), select=c(V1, V9, V10, V11, V12, V13, V14, V15))[,
  unique(V1), keyby=list(V9,V10,V11,V12,V13,V14,V15)
][, transmission_setting := seq(10,90,by=20)[V9+1] ][,
  scenario := scn_render(paste0(V10, V11, V12, V13, V14, V15))
][,
  list(seed=unique(V1)), by=list(transmission_setting, scenario)
], seed)

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

saveRDS(rbind(reduceditems, seroresults), "../longini2.RData")

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