#!/usr/bin/env Rscript
## calculate age distribution at time of vaccine; must be run in aggregated directory

#tar <- "~/Downloads/who-feb-2016/who-feb-2016-aggregated/"
#setwd(tar)

rm(list=ls(all.names = T))

args <- commandArgs(trailingOnly = T)
dbpath <- ifelse(is.na(args[1]), "~/Dropbox", args[1])

cat("dropbox path\n", dbpath,"\n")

bgfiles <- list.files(pattern = "^[0-4](_0)+_9(_0)+\\.")

if(length(bgfiles) == 0) stop("no source files in ", getwd())

translate_scenario <- function(dt) { # assumes transmission already separated
  ky <- key(dt)
  dt[grepl("_0$", scenario), scen := "coverageAT50%"]
  dt[grepl("^(0|1)_1", scenario), scen := "altVaccine" ]
  dt[grepl("^1_", scenario), scen := "catchUp"]
  dt[grepl("1_(0|1)$", scenario), scen := paste0(scen,"To30")]
  dt[grepl("0_0_1_9_0_1",scenario), scen := "reference"]
  dt[grepl("((0|1)_){3}[^9]+(_(0|1)){2}",scenario), scen := gsub("(_(1|0))+$","",gsub("^((1|0)_)+", "routineAT", scenario))]
  dt[grepl("^((0|1)_){2}0", scenario), scen := "noVaccine"]
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

require(reshape2)
require(data.table)

ref <- dcast.data.table(melt(
  rbindlist(lapply(bgfiles, fread, colClasses = c(scenario="character"))),
  id.vars = c("serial","scenario","age","outcome"), variable.name = "year", value.name = "count"
)[outcome != 0][,
  year := as.integer(gsub("y", "", year))
][,
  particle_id := floor(serial/800)
][,
  outcome := factor(c("mild","severe")[outcome], levels=c("mild","severe"), ordered = T)
], particle_id + scenario + age + year ~ outcome, value.var = "count")

ref[,
  symptomatic := mild + severe  
][,
  hospitalised := mild*phc + severe  
][,
  transmission_setting := seq(10,90,20)[as.integer(substr(scenario, 1, 1))+1]
][,
  scenario := substr(scenario,2,7)
]

translate_scenario(ref)

tmp <- ref[, {
  list(age=age, hospitalised_cases=cumsum(hospitalised), symptomatic_cases=cumsum(symptomatic))
}, keyby=list(scenario, particle_id, transmission_setting, year)]

#tmp[age == 100 & symptomatic_cases == 0, symptomatic_cases := 1]
#tmp[age == 100 & hospitalised_cases == 0, hospitalised_cases := 1]

mlt <- melt(tmp[,
  list(
    age=age,
    hospitalised_cases=hospitalised_cases/hospitalised_cases[.N],
    symptomatic_cases=symptomatic_cases/symptomatic_cases[.N]
  ),
  keyby=list(scenario, particle_id, transmission_setting, year)
][!is.na(hospitalised_cases) & !is.na(symptomatic_cases)], measure.vars = c("hospitalised_cases", "symptomatic_cases"))

combo <- mlt[,{
    mn = mean(value)
    se = sd(value)/sqrt(.N)
    list(value=mn, CI_low=mn-se, CI_high=mn+se)
  },
  keyby=list(scenario, transmission_setting, age, variable)  
][, outcome := gsub("_"," ",variable)][, group := "UF"][, outcome_denominator := "cumulative_proportion_at_baseline" ][, year:=0 ]

scens <- c("catchup","coverageAT50","reference","routineAT16")
fin <- rbindlist(lapply(scens, function(scen) copy(combo)[, scenario := scen ]))

saveRDS(fin, paste0(dbpath,"/CMDVI/Phase II analysis/Data/UF-Longini/longini-baseline-age-distro.rds"))

stop("run code in rstudio if you want to see plots")

## translate into symptomatic + hospitalised (overlapping)

require("ggplot2")
require("tidyr")
require("dplyr")
require("scales")

cbPalette <- c("Hopkins/UF"="#999999",
               "Imperial"="#E69F00",
               "Duke"="#56B4E9",
               "Notre Dame"="#009E73",
               "UF"="#F0E442",
               "Exeter/Oxford"="#0072B2",
               "Sanofi Pasteur"="#D55E00",
               "UWA"="#CC79A7",
               "A"="white",
               "Longini-ODE"="black"); 

noGroup_cbPalette <- c("#999999","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","black")

# cumulatage distribution of cases at baseline (assume that value is for age at birthday
# e.g.: age=4, value=0.5 -> 50% of cases are <4 yrs old) 
df_tmp=subset(fin, outcome %in% c("symptomatic cases","hospitalised cases"))
df_tmp$age=as.numeric(as.character(df_tmp$age))

ggplot(df_tmp, aes(x=age, y=value, group=scenario, color=scenario)) +
  geom_step(direction="hv") +
  facet_grid(outcome~transmission_setting, scale="free_y") +
  xlab("Age (yrs)") +
  ylab("cumulative proportion at baseline") +
  theme_bw() +
  scale_color_manual(values=noGroup_cbPalette) +
  coord_cartesian( xlim =c(0,60))

# age up to which 50% of cases are reported (not quite accurate for models with age groups)
tmp <- df_tmp %>% group_by(group, transmission_setting, scenario, outcome)
Age50percentCases <- tmp %>% summarise(minAge=min(age[value>0.5])-1) 
ggplot(subset(Age50percentCases, minAge<100),aes(x=group, y=minAge, color=scenario, fill=scenario)) +
  geom_bar(stat="identity", alpha=0.5, color=NA, position = "dodge") +
  facet_grid(outcome~transmission_setting, scale="free_y") +
  xlab("") +
  ylab("Age up tp which 50% of cases occur") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  scale_fill_manual(values=noGroup_cbPalette) + scale_color_manual(values=noGroup_cbPalette) 
