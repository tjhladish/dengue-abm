## economics calcs

tar <- "~/Downloads/who-feb-2016/who-feb-2016-aggregated/"
setwd(tar)

rm(list=ls(all.names = T))

require(reshape2)
require(data.table)
require(parallel)

poppath <- "/Volumes/Data/workspaces/dengue/pop-merida/pop-merida/population-merida.txt" # needs to be merida instead
age_cats <- c("<9yrs","9-18yrs","19+yrs","overall")

rr_severity_sec_vs_pri = 20
ph1 = 0.111
ph2 = ph1*1.88
phs = 1
ps1 = ph1*(1-1.88) / (-(rr_severity_sec_vs_pri-1) - ph1*1.88 + rr_severity_sec_vs_pri*ph1)
# ps2 = 20*ps1 = 0.1149942
# pc1 = 1 - ps1
# pc2 = 1 - ps2
phc = (0.111 - ps1)/(1-ps1)
cfr <- .00078

pop <- fread(poppath)[,list(pop=.N),keyby=age][,
  age_category := factor(ifelse(
    age < 9, age_cats[1], ifelse(
    age < 19, age_cats[2],
    age_cats[3]
  )), levels = age_cats, ordered = TRUE)
]

srcfiles <- list.files(pattern = "^[0-4].+\\.")

econdata <- subset(merge(
  dcast.data.table(melt(
    rbindlist(lapply(srcfiles, fread, colClasses = c(scenario="character"))),
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

# econfinal <- subset(econdata, select=c(particle_id, scen, transmission_setting, year, age, vac, amb, hosp, death))
# setkey(econfinal, particle_id, scen, transmission_setting, year, age)

econdata <- subset(econdata, select=-c(mild, severe, symptomatic, scenario, serial))

load("~/Downloads/lifeTablesMX.Rdata")
## from https://en.wikipedia.org/wiki/Demographics_of_Mexico 28 Jan 2016
# At birth	1.04 male(s)/female
# Under 15	1.05 male(s)/female
# 15–64 years	0.94 male(s)/female
# 65 and over	0.81 male(s)/female
genderratio <- data.table(
  ms=c(1.04, rep(1.05, 14), rep(0.94, 64-15+1), rep(0.81, 100-65+1)))[,
  age:=0:(.N-1)
][,
  `:=`(
    mp = ms/(ms+1),
    fp = 1/(ms+1)
  )
][,list(mp, fp), keyby=age]

age.max <- 100
mort <- data.table(femaleLT)[,prob.die,keyby=list(age=min.age)][data.table(maleLT)[,prob.die,keyby=list(age=min.age)]][genderratio][,list(p.death=fp*prob.die+mp*i.prob.die),keyby=age]$p.death
mort[age.max+1] <- 1

#Based on...but with substantial refactoring
#Cost-effectiveness calculations for WHO dengue vaccine modelling consortium
#(c) Mark Jit and Stefan Flasche 2016
#mark.jit@lshtm.ac.uk; January 2016; version 2
#changelog v2
#*Added societal perspective scenarios for reference case
#*Changed main outcome from ICER to net monetary benefit
#*Catchup now compared to reference rather than no vacc
#*Threshold for CE set to 2000, 4000, ..., 10000

require(tidyr)
require(dplyr)

#############################################
# input from modellers
#############################################
# n.group="UF"

# #############################################
# # Economic parameters and scenarios
# #############################################
# 
#   econ.scenarios=c("PHL","nodisc","society")

disc.cost = 0.03	#discount rate for costs
disc.daly = 0.03	#discount rate for dalys
persp = 0			#perspective 0=health care provider, 1=society
thresholds=(0:5)*2000  #threshold for cost-effectiveness

econ.para.bra=list(
  c=list(
    cvac = 30.50,	#HCP cost of vaccination
    camb = 60,		#HCP cost of treating VCD case
    chosp = 200,		#HCP cost of treating hospitalised VCD case
    cdeath= 0		#HCP cost of dengue death
  ),
  
  sc=list(
    cvac = 30.50,	#societal cost of vaccination
    camb = 200,		#societal cost of treating VCD case
    chosp = 500,	#societal cost of treating hospitalised VCD case
    cdeath = 11000 #societal cost of dengue death
  ),	
  
  d=list(
    d.vac = 0,			#dalys due to vaccination
    d.amb = 0.545*4/365,	#dalys due to ambulatory case
    d.hosp = 0.545*14/365	#dalys due to hospitalised case
  )
)

econ.para.phl=list(
  c=list(
    cvac = 30.50,	#HCP cost of vaccination
    camb = 20,		#HCP cost of treating VCD case
    chosp = 400,		#HCP cost of treating hospitalised VCD case
    cdeath = 0		#HCP cost of dengue death
  ),
  sc=list(
    cvac = 30.50,	#societal cost of vaccination
    camb = 40,		#societal cost of treating VCD case
    chosp = 500,	#societal cost of treating hospitalised VCD case
    cdeath = 3000	#societal cost of treating hospitalised VCD case
  ),
  d=list(
    d.vac = 0,		#dalys due to vaccination
    d.amb = 0.545*4/365,	#dalys due to ambulatory case
    d.hosp = 0.545*14/365	#dalys due to hospitalised case
  )
)

econ.scenarios=
  c("PHL","nodisc","society",
    "amb.lo"  , "amb.hi"  , "hosp.lo"  , "hosp.hi"  , "death.lo"  , "death.hi",
    "amb.lo.s", "amb.hi.s", "hosp.lo.s", "hosp.hi.s", "death.lo.s", "death.hi.s")

econ.para.bra.amb.lo=econ.para.bra 
econ.para.bra.amb.hi=econ.para.bra
econ.para.bra.hosp.lo=econ.para.bra 
econ.para.bra.hosp.hi=econ.para.bra 
econ.para.bra.death.lo=econ.para.bra
econ.para.bra.death.hi=econ.para.bra

econ.para.bra.amb.lo$c$camb=econ.para.bra$c$camb*0.5
econ.para.bra.amb.hi$c$camb=econ.para.bra$c$camb*1.5
econ.para.bra.hosp.lo$c$chosp=econ.para.bra$c$chosp*0.5
econ.para.bra.hosp.hi$c$chosp=econ.para.bra$c$chosp*1.5
econ.para.bra.death.lo$c$cdeath=econ.para.bra$c$cdeath*0.5
econ.para.bra.death.hi$c$cdeath=econ.para.bra$c$cdeath*1.5

econ.para.bra.amb.lo$sc$camb=econ.para.bra$sc$camb*0.5
econ.para.bra.amb.hi$sc$camb=econ.para.bra$sc$camb*1.5
econ.para.bra.hosp.lo$sc$chosp=econ.para.bra$sc$chosp*0.5
econ.para.bra.hosp.hi$sc$chosp=econ.para.bra$sc$chosp*1.5
econ.para.bra.death.lo$sc$cdeath=econ.para.bra$sc$cdeath*0.5
econ.para.bra.death.hi$sc$cdeath=econ.para.bra$sc$cdeath*1.5


econ.para.bra.amb.lo$d$d.amb=econ.para.bra$d$d.amb*0.5
econ.para.bra.amb.hi$d$d.amb=econ.para.bra$d$d.amb*1.5
econ.para.bra.hosp.lo$d$d.hosp=econ.para.bra$d$d.hosp*0.5
econ.para.bra.hosp.hi$d$d.hosp=econ.para.bra$d$d.hosp*1.5

############################################# 
# calcuate outcomes
############################################# 

calcLifExp = function(age.max, mort, disc.daly){
  surv=rep(0,age.max+1)		#survival
  disc.surv=rep(0,age.max+1)	#discounted survival
  life.exp=rep(0,age.max+1)	#life expectancy
  disc.life.exp=rep(0,age.max+1)#discounted life expectancy
  surv[1]=1
  disc.surv[1]=1
  for(a in 1:age.max){
    surv[a+1] = surv[a] * (1-mort[a])
    disc.surv[a+1] = disc.surv[a] * (1-mort[a])/(1+disc.daly)
  }
  for(a in 0:age.max){
    life.exp[a+1]=sum(surv[(a+1):(age.max+1)])/surv[a+1]
    disc.life.exp[a+1]=sum(disc.surv[(a+1):(age.max+1)])/disc.surv[a+1]
  }
  return(cbind(survival=disc.surv, lifeExp=disc.life.exp))
}

out=calcLifExp(age.max, mort, disc.daly)
#  disc.surv=out[,"survival"] # unused
disc.life.exp=out[,"lifeExp"]

year.max=econdata[,max(year)]

disc.cost.vec = 1/rep(disc.cost+1, year.max+1)^(0:year.max)
disc.daly.vec = 1/rep(disc.daly+1, year.max+1)^(0:year.max)
nodisc = rep(1, length(disc.cost.vec))

precomputeCost <- function(tardt, costsrc, wh="c") with(c(costsrc[[wh]], costsrc$d), {
  tardt[,`:=`(
    n.vac = vac,
    cost.vac = vac * cvac,
    cost.treat = (amb*camb + hosp*chosp + death*cdeath),
    daly.vac = vac * d.vac,
    daly.treat = (amb*d.amb + hosp*d.hosp),
    daly.death = death
  )]
})

slice <- function(base) base[,list(n.vac, cost.vac, cost.treat, daly.vac, daly.treat, daly.death), keyby=list(thresh, transmission_setting, particle_id, age, year, scen)]
extract_ref <- function(dt, wh) subset(dt[scen == wh], select=-scen)

# for each of these...
econ.scenarios=
  list(
    "BRA"=econ.para.bra,
    "PHL"=econ.para.phl,
    "nodisc"=econ.para.bra, "society"=econ.para.bra, # but with different other inputs
    "amb.lo"=econ.para.bra.amb.lo, "amb.hi"=econ.para.bra.amb.hi,
    "hosp.lo"=econ.para.bra.hosp.lo, "hosp.hi"=econ.para.bra.hosp.hi,
    "death.lo"=econ.para.bra.death.lo  , "death.hi"=econ.para.bra.death.hi,
    "amb.lo.s"=econ.para.bra.amb.lo, "amb.hi.s"=econ.para.bra.amb.hi,
    "hosp.lo.s"=econ.para.bra.hosp.lo, "hosp.hi.s"=econ.para.bra.hosp.hi,
    "death.lo.s"=econ.para.bra.death.lo  , "death.hi.s"=econ.para.bra.death.hi)

bases=list("c","c","c","sc","c","c","c","c","c","c","sc","sc","sc","sc","sc","sc")
regs=c(list("BRA","PHL"),rep("BRA",length(econ.scenarios)-2))
discounts=c(list(F,F,T),rep(T,length(econ.scenarios)-3))

# need one of these
econfinalBRA <- econdata[econdata[,list(thresh=thresholds),keyby=key(econdata)], allow.cart=T]
setkeyv(econfinalBRA, c("thresh", key(econdata)))

fun <- function(para, basis, reg, discounting, scen.name) {
  init <- copy(econfinalBRA)
  cont <- slice(precomputeCost(init, para, basis))
  refed <- extract_ref(econfinalBRA, "reference")
  novaced <- extract_ref(econfinalBRA, "noVaccine")
  joins(cont, refed, novaced, reg, ifelse(basis=="c","ind","soc"), scen.name,
        if(discounting) nodisc else disc.cost.vec, if(discounting) nodisc else disc.daly.vec,
        discounting
  )
}

calc <- function(.SD, cst.disc, d.disc, thresh) with(.SD, {
  net.cost=sum((cost.vac-i.cost.vac+cost.treat-i.cost.treat)*cst.disc[year+1])
  net.daly=sum((daly.vac-i.daly.vac+daly.treat-i.daly.treat+(daly.death-i.daly.death)*disc.life.exp[age+1])*d.disc[year+1])
  net.cost.treat=sum((cost.treat-i.cost.treat)*cst.disc[year+1])
  net.n.vac=sum((n.vac-i.n.vac)*cst.disc[year+1])
  icer = -net.cost/net.daly	#Incremental cost per DALY averted
  nmb = -(net.cost+net.daly*thresh) 
  nmbpp = nmb/net.n.vac
  threshold.cost = -(net.cost.treat + net.daly * thresh)/net.n.vac	#Threshold cost per person vaccinated
  
  # manually correct for rounding errors
  if(abs(net.cost)<.1) icer=0
  if(abs(net.cost)<.1 & abs(net.daly)<.1) threshold.cost=0
  
  list(Cost=net.cost, DALY=net.daly, ICER=icer, ThresholdCosts=threshold.cost, NMB=nmb, NMBpp=nmbpp)
})

quan <- function(var, nm, n) {
  mn = mean(var)
  se = 2*sd(var)/sqrt(n)
  qres = c(mn, mn-se, mn+se)
  names(qres) <- paste0(c("value","CI_low.","CI_high."), nm)
  as.list(qres)
}

agg <- function(base, ref, discc, discd, ...) base[...][ref][,
  calc(.SD, discc, discd, thresh),
  by=list(thresh, scen, transmission_setting, particle_id)
][,{
    qThresholds <- quan(ThresholdCosts, "thresh", .N)
    qNMBpp <- quan(NMBpp, "nmbpp", .N)
    qICER <- quan(NMBpp, "icer", .N)
    qNMB <- quan(NMB, "nmb", .N)
    c(
      list(value.thresh=mean(ThresholdCosts), value.nmbpp=mean(NMBpp), value.NMB=mean(NMB), value.icer=mean(ICER)),
      qThresholds, qNMBpp, qICER, qNMB
    )
  },
  by=list(thresh, scen, transmission_setting)
]

joins <- function(base, ref, noVac, region, typ, eco, cost.vec, daly.vec, discounting) {
  rbind(
    agg(base, ref, cost.vec, daly.vec, scen %in% c("catchUp","catchupTo30"))[, disc := discounting ],
    agg(base, noVac, cost.vec, daly.vec, !(scen %in% c("catchUp","catchupTo30","noVaccine")))[, disc := discounting ]
  )[,
    `:=`(reg = region, cost = typ, ecoscn = eco)
  ]
}

allres <- rbindlist(mcmapply(
  fun,
  para=econ.scenarios,
  basis=bases,
  reg=regs,
  discounting=discounts,
  scen.name = names(econ.scenarios),
  mc.cores = detectCores()-1
)
#   joins(econfinalBRAInd, econfinalBRAIndref, econfinalBRAIndnoVac, "BRA", "ind"),
#   joins(econfinalBRASoc, econfinalBRASocref, econfinalBRASocnoVac, "BRA", "soc"),
#   joins(econfinalPHLInd, econfinalPHLIndref, econfinalPHLIndnoVac, "PHL", "ind"),
#   joins(econfinalPHLSoc, econfinalPHLSocref, econfinalPHLSocnoVac, "PHL", "soc")
)[, year := paste0("cum", year.max+1) ][, age := "overall" ][, group := "UF" ][,
  scenario := paste(scen,reg,cost,ecoscen,sep="-")
]

# precomputeCost(econfinalBRAInd, econ.para.bra, "c")
# precomputeCost(econfinalBRASoc, econ.para.bra, "sc")
# precomputeCost(econfinalPHLInd, econ.para.phl, "c")
# precomputeCost(econfinalPHLSoc, econ.para.phl, "sc")
# 
# econfinalBRAInd <- slice(econfinalBRAInd)
# econfinalBRASoc <- slice(econfinalBRASoc)
# econfinalPHLInd <- slice(econfinalPHLInd)
# econfinalPHLSoc <- slice(econfinalPHLSoc)
# 
# econfinalBRAIndref <- extract_ref(econfinalBRAInd, "reference")
# econfinalBRAIndnoVac <- extract_ref(econfinalBRAInd, "noVaccine")
# econfinalBRASocref <- extract_ref(econfinalBRASoc, "reference")
# econfinalBRASocnoVac <- extract_ref(econfinalBRASoc, "noVaccine")
# 
# econfinalPHLIndref <- extract_ref(econfinalPHLInd, "reference")
# econfinalPHLIndnoVac <- extract_ref(econfinalPHLInd, "noVaccine")
# econfinalPHLSocref <- extract_ref(econfinalPHLSoc, "reference")
# econfinalPHLSocnoVac <- extract_ref(econfinalPHLSoc, "noVaccine")



econres <- rbind(
  allres[,
    list(outcome=paste0("NMBpp-",thresh), value=value.nmbpp, CI_low=CI_low.nmbpp, CI_high=CI_high.nmbpp, outcome_denominator=NA, year, age="overall"),
    by=list(transmission_setting, group, scenario)
  ],
  allres[,
    list(outcome=paste0("NMB-",thresh), value=value.nmbpp, CI_low=CI_low.nmbpp, CI_high=CI_high.nmbpp, outcome_denominator=NA, year, age="overall"),
    by=list(transmission_setting, group, scenario)
  ],
  allres[,
    list(outcome=paste0("ThreshholdCost-",thresh), value=value.thresh, CI_low=CI_low.thresh, CI_high=CI_high.thresh, outcome_denominator=NA, year, age="overall"),
    by=list(transmission_setting, group, scenario)
  ],
  allres[,
    list(outcome="ICER", value=unique(value.icer), CI_low=unique(CI_low.icer), CI_high=unique(CI_high.icer), outcome_denominator=NA, year=unique(year), age="overall"),
    by=list(transmission_setting, group, scenario)
  ]
)

setcolorder(econres, c("transmission_setting","scenario","year","outcome","age","value","CI_low","CI_high","outcome_denominator","group"))

saveRDS(econres,"~/Dropbox/CMDVI/Phase II analysis/Data/UF-Longini/longini-econ.rds")

stop("run script in rstudio if you want to pick out plots")

df <- within(econres,{
  scenario=gsub("-BRA-soc-disc","-society",scenario)
  scenario=gsub("-PHL-ind-disc","-PHL",scenario)
  scenario=gsub("-BRA-ind-nodisc","-nodisc",scenario)
  scenario=gsub("-BRA-ind-disc","",scenario)
  outcome = gsub("ThreshholdCost", "ThresholdPrice" ,outcome)
})
# 
# 
df_tmp=subset(df, (grepl("NMB-",outcome) | grepl("NMBpp-",outcome) | grepl("ThresholdPrice-",outcome)) & (year=="cum30") & (scenario == "reference"))
df_tmp$Threshold=as.numeric(matrix(unlist(strsplit(as.character(df_tmp$outcome),"-")),2)[2,])
df_tmp$outcome=(matrix(unlist(strsplit(as.character(df_tmp$outcome),"-")),2)[1,])
# 
# cbPalette <- c("Hopkins/UF"="#999999",
#                "Imperial"="#E69F00",
#                "Duke"="#56B4E9",
#                "Notre Dame"="#009E73",
#                "UF"="#F0E442",
#                "Exeter/Oxford"="#0072B2",
#                "Sanofi Pasteur"="#D55E00",
#                "UWA"="#CC79A7",
#                "A"="white"); 
# 
noGroup_cbPalette <- c("#999999","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","black")
# 
ggplot(df_tmp,aes(x=Threshold, y=value, ymin=CI_low, ymax=CI_high, fill=scenario, color=scenario, group=scenario))+
  geom_point(alpha=0.6) +
  geom_line(alpha=0.6) +
  geom_ribbon(alpha=0.2, color=NA) +
  facet_grid(outcome~transmission_setting, scale="free_y") +
  xlab("Cost-effectiveness threshold (Costs per DALY)") +
  ylab("Net monetary benefit") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  scale_color_manual(values=noGroup_cbPalette)+ scale_fill_manual(values=noGroup_cbPalette) 

df_tmp=subset(df, (scenario %in% c("reference","reference-PHL", "reference-nodisc","reference-society")) & 
                (age=="overall") & 
                (year=="cum30") &
                (outcome %in% c("ThresholdPrice-2000"))) 

ggplot(df_tmp, aes(x=group, y=value, ymin=CI_low, ymax=CI_high, fill=scenario, color=scenario)) +
  geom_pointrange(alpha=0.6,position=position_dodge(.5)) +
  facet_grid(outcome~transmission_setting, scale="free_y") +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  scale_color_manual(values=noGroup_cbPalette) 
