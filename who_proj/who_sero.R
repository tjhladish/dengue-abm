## seroprev data for who

#seropath <- "~/Downloads/who-feb-2016/processed_logs/" # path to serostatus processed logs - run process
#setwd(seropath)

rm(list=ls())

args <- commandArgs(trailingOnly = T)
poppath <- ifelse(is.na(args[1]), "~/git/dengue/pop-merida/pop-merida/population-merida.txt", args[1])
dbpath <- ifelse(is.na(args[2]), "~/Dropbox", args[2])

require(data.table)
require(reshape2)

cat("loading pop from: ",poppath,"\n")
cat("dropbox path: ", dbpath,"\n")

# age_cats <- c("<9yrs","9-18yrs","19+yrs","overall")
# 
# pop <- fread(poppath)[,list(pop=.N),keyby=age][,
#   age_category := factor(ifelse(
#     age < 9, age_cats[1], ifelse(
#      age < 19, age_cats[2],
#      age_cats[3]
#     )), levels = age_cats, ordered = TRUE
#   )
# ]

seroscn <- setkey(subset(fread("seroscenarios.csv"), select=c(V1, V9, V10, V11, V12, V13, V14, V15))[,
  unique(V1), keyby=list(V9,V10,V11,V12,V13,V14,V15)
][, transmission_setting := seq(10,90,by=20)[V9+1] ][,
  scenario := paste0(V10, V11, V12, V13, V14, V15)
][,
  list(seed=unique(V1)), keyby=list(transmission_setting, scenario)
], seed)

translate_scenario <- function(dt) { # assumes transmission already separated
  ky <- key(dt)
  dt[substr(scenario,6,6) == "0", scen := "coverageAT50"]
  dt[substr(scenario,2,2) == "1", scen := "altVaccine" ]
  dt[substr(scenario,1,1) == "1", scen := "catchUp"]
  dt[substr(scenario,5,5) == "1", scen := paste0(scen,"To30")]
  dt[grepl("^0{2}10{2}1",scenario), scen := "reference"]
  dt[substr(scenario,4,4) == "1", scen := "routineAT16"]
  dt[substr(scenario,3,3) == "0", scen := "noVaccine"]
  setkeyv(subset(dt[, scenario:=scen ], select=-scen), ky)
}

translate_scenario(seroscn)

serosrv <- setkey(fread("seroprevalence.csv")[,
  list(seed=unique(seed)), keyby=list(year, vax_age, seropositive, pop)
][, year := year-50 ], seed, year)

preseroresults <- setnames(subset(serosrv[seroscn], select=-c(scen,seed)), "vax_age", "age")
setkey(preseroresults, age, scenario, transmission_setting, year)

seroresults <- preseroresults[,{
    mn <- mean(seropositive/pop)
    se <- sqrt(mn*(1-mn)/.N)
    list(value=mn, CI_low=mn-se, CI_high=mn+se)
  },
  keyby=list(scenario, transmission_setting, year, age)
]

seroresults[, age := paste0(age, "yo") ][, outcome := "seropositive" ][, group:="UF" ][, outcome_denominator:="proportion" ]

saveRDS(seroresults[age=="9yo"],paste0(dbpath, "/CMDVI/Phase II analysis/Data/UF-Longini/longini-serostatus.rds"))

stop()

require(ggplot2)

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
scenarios=c("reference","catchUp","routineAT16","coverageAT50","noVaccine")
for (scen in scenarios){
  df=subset(seroresults, (scenario==scen) & (age=="9yo") )
  df_tmp=subset(df, outcome=="seropositive")
  df_tmp$year=as.numeric((as.character(df_tmp$year)))
  
  p<-ggplot(df_tmp, aes(x= year, y=value, ymin=CI_low, ymax=CI_high, fill=group, color=group, group=group)) +
    geom_line() + ggtitle(scen) +
    geom_ribbon(alpha=0.2, color=NA) +
    facet_grid(.~transmission_setting, scale="free_y") +
    xlab("time after introduction of CYD (years)") +
    ylab("proportion seropositive\nat 9yrs of age") +
    theme_bw() + scale_y_continuous(breaks=c(0,10,30,50,70,90,100)/100, limits=c(0,1)) +
    scale_fill_manual(values=cbPalette) + scale_color_manual(values=cbPalette)
  print(p)
}