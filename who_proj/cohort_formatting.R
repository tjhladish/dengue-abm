## cohort review

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

tar <- "~/Downloads/who-feb-2016/who-feb-2016-aggregated/"
setwd(tar)

vacfiles <- list.files(pattern = "^[0-4][01]{2}1[01]{3}\\.")
bgfiles <- list.files(pattern = "^[0-4]0+\\.")

seedkeyFG <- setkey(rbindlist(lapply(vacfiles, function(nm) {
  c(list(seed=as.numeric(gsub(".+\\.","", nm))), fread(nm, nrows = 1, colClasses = c(scenario="character"))[,list(serial, scenario)])
})), seed)[, particle_id := floor(serial/80) ]

seedkeyBG <- setkey(rbindlist(lapply(bgfiles, function(nm) {
  c(list(seed=as.numeric(gsub(".+\\.","", nm))), fread(nm, nrows = 1, colClasses = c(scenario="character"))[,list(serial, scenario)])
})), seed)[, particle_id := floor(serial/80) ]

readInCohort <- function(src, sk, cnames = c("seed","year","serostatus","vaccination","severity","count")) setkey(subset(
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

cohortdata <- readInCohort("../cohort/cohort.vaccinations", seedkeyFG)
refcohortdata <- readInCohort("../cohort/cohort.no-vaccinations", seedkeyBG,c("seed","year","serostatus","vaccination","severity","count","routineage"))

process_severity <- function(tar) {
  tar[severity!=0, infections := count]
  tar[severity>1, symptomatic_cases := count ]
  tar[severity>2, hospitalised_cases := count ]
  tar[severity==2, c("hospitalised_cases","deaths") := {
    h = rbinom(.N, count, phc)
    d = rbinom(.N, h, cfr)
    list(h,d)
  }]
  tar[severity==3, deaths := rbinom(.N, hospitalised_cases, cfr)]
}

cohortdata <- translate_scenario(process_severity(cohortdata))
refcohortdata <- translate_scenario(process_severity(refcohortdata))