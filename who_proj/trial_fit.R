#!/usr/bin/env Rscript
args <- c("~/Downloads/who_revision_data2/who-80-rev-forecast-attempt2.sqlite")
## read in sim outputs
## out format: Groups Trial    Description    Arm       Strata        mean       lower       upper

rm(list=ls())

args <- commandArgs(trailingOnly = T)
sqlpath <- ifelse(is.na(args[1]), "~/Downloads/cyd14_match-correct_mechanism-extra_year.sqlite", args[1])
store <- !is.na(args[2])
dbpath <- ifelse(is.na(args[3]), "~/Dropbox", args[3])

cat("loading from: ",sqlpath,"\n")
cat("storing? ", store, "\n")
cat("dropbox path\n", dbpath,"\n")

require(RSQLite)
require(data.table)
require(ggplot2)

getMetrics <- function(sqlite) {
  db <- RSQLite::dbConnect(RSQLite::SQLite(), sqlite)
  result <- data.table(RSQLite::dbGetQuery(db,
    "select m.*, p.seed from metrics m, parameters p, jobs j where m.serial = p.serial and j.serial = m.serial and j.status == 'D';"
  ))
  RSQLite::dbDisconnect(db)
  result
}

#init <- subset(getMetrics(sqlpath)[beta_multiplier == 4], select=-beta_multiplier)[,CYD:=14]
init <- getMetrics(sqlpath)

ps1 = 0.005749711
phc = (0.111 - ps1)/(1-ps1)

init[,
  `:=`(
    hosp_vaccine_12_14 = hosp_severe_vaccine_12_14+ hosp_mild_vaccine_12_14*phc,
    hosp_vaccine_2_5   = hosp_severe_vaccine_2_5  + hosp_mild_vaccine_2_5*phc,
    hosp_vaccine_6_11  = hosp_severe_vaccine_6_11 + hosp_mild_vaccine_6_11*phc,
    hosp_placebo_12_14 = hosp_severe_placebo_12_14+hosp_mild_placebo_12_14*phc,
    hosp_placebo_2_5   = hosp_severe_placebo_2_5  +hosp_mild_placebo_2_5*phc,
    hosp_placebo_6_11  = hosp_severe_placebo_6_11 +hosp_mild_placebo_6_11*phc,
    hosp_vaccine_12_14_y2 = hosp_severe_vaccine_12_14_y2+ hosp_mild_vaccine_12_14_y2*phc,
    hosp_vaccine_2_5_y2   = hosp_severe_vaccine_2_5_y2  + hosp_mild_vaccine_2_5_y2*phc,
    hosp_vaccine_6_11_y2  = hosp_severe_vaccine_6_11_y2 + hosp_mild_vaccine_6_11_y2*phc,
    hosp_placebo_12_14_y2 = hosp_severe_placebo_12_14_y2+ hosp_mild_placebo_12_14_y2*phc,
    hosp_placebo_2_5_y2   = hosp_severe_placebo_2_5_y2  + hosp_mild_placebo_2_5_y2*phc,
    hosp_placebo_6_11_y2  = hosp_severe_placebo_6_11_y2 + hosp_mild_placebo_6_11_y2*phc
  )
]

init <- subset(init, select=grep("(severe|mild)", names(init), value = T, invert = T))
init[,CYD:=14]

mlt <- subset(melt(init, id.var=c("seed","serial","CYD")), select=-seed)
# need to create: y1+y2 result
y1dat <- mlt[grepl("hosp", variable) & !grepl("y2", variable)]
y2dat <- mlt[grepl("y2", variable)]
y2dat[, variable := gsub("_y2","",variable)]
y12dat <- merge(y1dat, y2dat, by=c("serial","CYD","variable"))
y12dat[, value:=value.x+value.y]
y12dat$value.x <- y12dat$value.y <- NULL
y12dat[, variable := gsub("hosp","hosp12",variable)]
rrext <- function(dat,newnm) {
  pl <- dat[grepl("placebo", variable)]
  pl[, variable := gsub("(placebo|vaccine)_","",variable)]
  va <- dat[grepl("vaccine", variable)]
  va[, variable := gsub("(placebo|vaccine)_","",variable)]
  RR <- merge(va, pl, by=c("serial","CYD","variable"))
  RR[, list(
    variable=gsub("hosp",newnm,variable),
    value=ifelse(value.x == 0, 0, value.x/value.y)
  ), by=list(serial, CYD)]
}
# y1 vaccine / placebo
# y2 vaccine / placebo
y1RR <- rrext(y1dat,"hosp1RR")
y2RR <- rrext(y2dat,"hosp2RR")

mlt <- rbind(
  mlt[!grepl("y2", variable)],
  y12dat,
  y1RR[!is.infinite(value)],
  y2RR[!is.infinite(value)]
)

wider <- mlt[,{
  mn = mean(value)
  se = 2*sd(value)/sqrt(.N)
  list(mean=mn, lower=mn-se, upper=mn+se, med=median(value),
       qlo=quantile(value,probs = 0.25), qlo=quantile(value,probs = 0.75))
}, keyby=list(Trial = paste0("CYD", CYD), variable)][!grepl("infec", variable)]

wider[, Arm := "Seronegative"]
wider[grepl("placebo",variable), Arm:="Placebo"]
wider[grepl("vaccine",variable), Arm:="Vaccine"]
wider[grepl("RR",variable), Arm:="RR"]
wider[, Strata := gsub(".+_(\\d+)_(\\d+)$","\\1-\\2y", variable) ]

wider[grepl("(neg|pos)",Strata),
  Strata := ifelse(grepl("neg",Strata),"Seronegative","Seropositive")
]
wider[, Description := "Full_trial"]
wider[grepl("hosp_", variable), Description := "LTFUY1"]
wider[grepl("hosp12_", variable), Description := "LTFUY12"]
wider[grepl("hosp1RR_", variable), Description := "LTFUY1_RR"]
wider[grepl("hosp2RR_", variable), Description := "LTFUY2_RR"]
wider[grepl("Sero",Strata), Description := "Immuno"]
wider[grepl("Sero",Arm), Description := "Immuno_seropos"]

output <- subset(wider,select=-variable)[,Groups:="UF"][, N:=NA ][, Cases:=NA ]

outpath <- paste0(dbpath,"/CMDVI/Phase II analysis/Data/UF-Longini/longini-fit.rds")

if (store) saveRDS(output, outpath)

#create.plots(df15,"CYD15",agegps15)
df_data=read.csv(paste0(dbpath, "/CMDVI/Phase II analysis/Data/_Fit/FitData.Rdata.csv"))
df_UF <- output[,names(df_data), with=F]
df = rbind(df_data, df_UF)

cbPalette <- c("Hopkins/UF"="#999999",
               "Imperial"="#E69F00",
               "Duke"="#56B4E9",
               "Notre Dame"="#009E73",
               "UF"="#F0E442",
               "Exceter/Oxford"="#0072B2",
               "Sanofi Pasteur"="#D55E00",
               "UWA"="#CC79A7",
               "Data"="black"); 

agegps14=c("2-5y", "6-11y", "12-14y")
agegps15=c("9-11y", "12-16y")

create.plots<-function(df,name,agegps,save=T){
  
  # Proportion Seronegative
  df_tmp=subset(df, Description=="Immuno_seropos" & Strata!="All")
  df_tmp$Strata = factor(df_tmp$Strata,agegps)
  p=ggplot(df_tmp, aes(x=Strata, y=mean, ymin=lower, ymax=upper, color=Groups)) +
    geom_pointrange(position=position_dodge(.2)) +
    scale_y_continuous("proportion seronegative \nat vaccination") +
    expand_limits(y = 0) +
    scale_x_discrete("Age group") +
    theme_bw() + scale_color_manual(values=cbPalette) 
  if (save) {
    ggsave(paste(path_figures,name,"Proportion_Seronegative",plot_type,sep=""),p,dpi=dpi_out,units="cm",width=20, height=8)
  } else print(p)
  # Attack rate by serostatus
  df_tmp=subset(df, Description=="Immuno")
  p=ggplot(df_tmp, aes(x=Arm, y=mean, ymin=lower, ymax=upper, color=Groups)) +
    geom_pointrange(position=position_dodge(.3)) +
    facet_grid(.~Strata) +
    scale_y_continuous( "Attack rate") +
    expand_limits(y = 0) +
    scale_x_discrete("") +
    theme_bw() + scale_color_manual(values=cbPalette)
  if (save) {
    ggsave(paste(path_figures,name,"Attack_rate_sero",plot_type,sep=""),p,dpi=dpi_out,units="cm",width=20, height=8)
  } else print(p)
  # Attack rate by age
  df_tmp=subset(df, Description=="Full_trial")
  df_tmp$Strata = factor(df_tmp$Strata,agegps)
  p=ggplot(df_tmp, aes(x=Arm, y=mean, ymin=lower, ymax=upper, color=Groups)) +
    geom_pointrange(position=position_dodge(.3)) +
    facet_grid(.~Strata) +
    scale_y_continuous("Attack rate") +
    expand_limits(y = 0) +
    scale_x_discrete("") +
    theme_bw() + scale_color_manual(values=cbPalette)
  if (save) {
    ggsave(paste(path_figures,name,"Attack_rate_age",plot_type,sep=""),p,dpi=dpi_out,units="cm",width=20, height=8)
  } else print(p)
  
  
  # Hospital phase attack rate
  df_tmp=subset(df, Description=="LTFUY1")
  df_tmp$Strata = factor(df_tmp$Strata,agegps)
  p=ggplot(df_tmp, aes(x=Arm, y=mean, ymin=lower, ymax=upper, color=Groups)) +
    geom_pointrange(position=position_dodge(.2)) +
    facet_grid(.~Strata) +
    scale_y_continuous("Attack rate \n in passive hospital surveillance") +
    expand_limits(y = 0) +
    scale_x_discrete("") +
    theme_bw() + scale_color_manual(values=cbPalette)
  if (save) {
    ggsave(paste(path_figures,name,"Hospital_phase_attack_rate",plot_type,sep=""),p,dpi=dpi_out,units="cm",width=20, height=8)
  } else print(p)
} 

create.plots(df,"CYD14",agegps14, save = F)
