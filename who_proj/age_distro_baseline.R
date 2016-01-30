## calculate age distribution at time of vaccine

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

require(reshape2)
require(data.table)

tar <- "~/Downloads/who-feb-2016/who-feb-2016-aggregated/"
setwd(tar)

poppath <- "~/git/dengue/pop-merida/pop-merida/population-merida.txt"

bgfiles <- list.files(pattern = "^[0-4]0+\\.")

ref <- dcast.data.table(melt(
  rbindlist(lapply(bgfiles, fread, colClasses = c(scenario="character"))),
  id.vars = c("serial","scenario","age","outcome"), variable.name = "year", value.name = "count"
)[outcome != 0][,
  year := as.integer(gsub("y", "", year))
][year==0][,
  particle_id := floor(serial/80)
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
}, keyby=list(scenario, particle_id, transmission_setting)]

exclude <- setkey(tmp[age == 100 & symptomatic_cases == 0, list(particle_id, transmission_setting)], particle_id, transmission_setting)

tmp[!exclude]

# mlt <- melt([,
#   list(age=age, hospitalised_cases=hospitalised_cases/hospitalised_cases[.N], symptomatic_cases=symptomatic_cases/symptomatic_cases[.N]),
#   keyby=list(scenario, particle_id, transmission_setting)
# ], measure.vars = c("hospitalised_cases", "symptomatic_cases"))

mlt[,{
    mn = mean(value)
    se = sd(value)/sqrt(.N)
    browser()
    list(value=mn, CI_low=mn-se, CI_high=mn+se)
  },
  keyby=list(scenario, transmission_setting, age, variable)  
]

## translate into symptomatic + hospitalised (overlapping)

df_tmp=subset(df, outcome %in% c("symptomatic cases","hospitalised cases"))
df_tmp$age=as.numeric(as.character(df_tmp$age))
table(df_tmp$age>100,df_tmp$group)

p<-ggplot(df_tmp, aes(x=age, y=value, group=group, color=group)) +
  geom_step(direction="hv") +
  facet_grid(outcome~transmission_setting, scale="free_y") +
  xlab("Age (yrs)") +
  ylab("cumulative proportion at baseline") +
  theme_bw() +
  scale_color_manual(values=cbPalette) 