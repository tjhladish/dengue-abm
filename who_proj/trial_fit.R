## read in sim outputs
## out format: Groups Trial    Description    Arm       Strata        mean       lower       upper

rm(list=ls())

df_data=read.csv("~/Dropbox/CMDVI/Phase II analysis/Data/_Fit/FitData.Rdata.csv")

require(RSQLite)
require(data.table)

getMetrics <- function(sqlite) {
  db <- RSQLite::dbConnect(RSQLite::SQLite(), sqlite)
  result <- data.table(RSQLite::dbGetQuery(db, "select seed, * from metrics m, parameters p where m.serial = p.serial;"))
  RSQLite::dbDisconnect(db)
  result
}

init <- subset(getMetrics("~/Downloads/who-feb-2016/cyd14_match.sqlite")[beta_multiplier == 4], select=-c(serial, seed, seed, mild_EF, severe_EF, base_path, sec_sev, pss_ratio, exp_coef, num_mos, beta, beta_multiplier))

ps1 = 0.005749711
phc = (0.111 - ps1)/(1-ps1)

init[,
  hosp_vaccine_12_14 := hosp_severe_vaccine_12_14+hosp_mild_vaccine_12_14*phc
][,
  hosp_vaccine_2_5 := hosp_severe_vaccine_2_5+hosp_mild_vaccine_2_5*phc
][,
  hosp_vaccine_6_11 := hosp_severe_vaccine_6_11+hosp_mild_vaccine_6_11*phc
][,
  hosp_placebo_12_14 := hosp_severe_placebo_12_14+hosp_mild_placebo_12_14*phc
][,
  hosp_placebo_2_5 := hosp_severe_placebo_2_5+hosp_mild_placebo_2_5*phc
][,
  hosp_placebo_6_11 := hosp_severe_placebo_6_11+hosp_mild_placebo_6_11*phc
]

init <- subset(init, select=grep("(severe|mild)", names(init), value = T, invert = T))

mlt <- subset(melt(init, id.var=c("seed","serial","CYD")), select=-seed)

mlt[grepl("hosp", variable), value, by=list(serial, severity=gsub(".+(severe|mild).+","\\1",variable))]

wider <- mlt[,{
  mn = mean(value)
  se = 2*sd(value)/sqrt(.N)
  list(mean=mn, lower=mn-se, upper=mn+se, Groups="longini")
}, keyby=list(Trial = paste0("CYD", CYD), variable)]

wider[, Arm := "Seronegative"]
wider[grepl("placebo",variable), Arm:="Placebo"]
wider[grepl("vaccine",variable), Arm:="Vaccine"]
wider[, Strata := gsub("[^1]+(\\d+)_(\\d+)$","\\1-\\2y", variable) ]

wider[grepl("(neg|pos)",Strata),
  Strata := ifelse(grepl("neg",Strata),"Seronegative","Seropositive")
]
wider[, Description := "Full_trial"]
wider[grepl("hosp", variable), Description := "LTFUY1"]
wider[grepl("Sero",Strata), Description := "Immuno"]
wider[grepl("Sero",Arm), Description := "Immuno_seropos"]

saveRDS(subset(wider,select=-variable)[,Groups:="UF"][, N:=NA ][, Cases:=NA ], "~/Dropbox/CMDVI/Phase II analysis/Data/UF-Longini/longini-fit.rds")

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

#create.plots(df15,"CYD15",agegps15)
df_UF <- readRDS("~/Dropbox/CMDVI/Phase II analysis/Data/UF-Longini/longini-fit.rds")[,names(df_data), with=F]
df = rbind(df_data, df_UF)

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