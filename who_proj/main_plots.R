# written by Stefan Flasche


####################################
#user input
####################################

#path_code="D:/GitHub/CMDVI---Analysis/"
#path_code="C:\\Users\\Stefan Flasche\\Documents\\GitHub\\CMDVI---Analysis"
#path_data="Data/"
path_figures="../"
output <- "png"
#setwd("H:\\My Documents\\_Projects\\WHO\\WHO - Dengue\\AnalysisII\\")
#setwd("C:\\Users\\Stefan Flasche\\Documents\\Work\\WHO - Dengue\\AnalysisII")

####################################
# read libraries
####################################
require("ggplot2")
require("tidyr")
require("dplyr")
require("scales")

####################################
# load data
####################################
#load(paste(path_data,"data_mock_groupA.Rdata",sep="")); df_A=df
#load(paste(path_data,"data_mock_groupB.Rdata",sep="")); df_B=df
#load(paste(path_data,"data_mock_groupC.Rdata",sep="")); df_C=df
#df=rbind(df_A,df_B,df_C)

df <- readRDS("../longini2.RData")

####################################
# plots
####################################
  #colors
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  scn = "baseline"
  
  #vaccine impact
  for (my_year in c("cum10","cum30")){
    for (my_outcome_denominator in c("proportion averted","per 100,000 pop at risk")){
      
      df_tmp=subset(df, (age %in% c("overall", "<9yrs", "9-18yrs", "19+yrs")) & year==my_year &  outcome_denominator==my_outcome_denominator & (scenario == scn))
      dodge=position_dodge(width=15)
      df_tmp$outcome = gsub("_averted", "\naverted",df_tmp$outcome)
      df_tmp$outcome = gsub("_", " ",df_tmp$outcome)
 
      
      p=ggplot(df_tmp, aes(x=transmission_setting, y=value, ymin=CI_low, ymax=CI_high, fill=group, color=group)) +
        geom_bar(stat="identity", position=dodge, alpha=0.5) +
        geom_errorbar(stat="identity", position=dodge, alpha=0.8, width=0.25) +
        facet_grid(outcome~age, scale="free_y") +
        xlab("proportion of seropositive 9 years olds (transmission intensity)") +
        ylab(my_outcome_denominator) +
        theme_bw() +
        scale_fill_manual(values=cbPalette) + scale_color_manual(values=cbPalette) 
        
      #print(p)
      ggsave(paste(path_figures,"VaccineImpact_",my_outcome_denominator,"_",my_year,".",output,sep=""),p,dpi=300,units="cm",width=20, height=15)
    }
  }
    

  #vaccine effects in vaccinated cohort
  my_age= "cohort1"
  df_tmp=subset(df, age == my_age & (outcome_denominator=="per 100,000 pop at risk") & (scenario == scn))
  df_tmp$outcome = gsub("_averted", "\naverted",df_tmp$outcome)
  df_tmp$outcome = gsub("_", " ",df_tmp$outcome)
  df_tmp$year=as.numeric(as.character(df_tmp$year))
  
  p=ggplot(df_tmp, aes(x= year, y=value, ymin=CI_low, ymax=CI_high, fill=group, color=group, group=group)) +
    geom_line() +
    geom_ribbon(alpha=0.2, color=NA) +
    facet_grid(outcome~transmission_setting, scale="free_y") +
    xlab("time after introduction of CYD (years)") +
    ylab("per 100,000 pop at risk") +
    theme_bw() +
    scale_fill_manual(values=cbPalette) + scale_color_manual(values=cbPalette) 
  #print(p)
  ggsave(paste(path_figures,"VaccineEffectsInVaccCohort_",my_age,".",output,sep=""),p,dpi=300,units="cm",width=20, height=14)

  
  #change in vaccine impact over time
  df_tmp=subset(df, (age == "overall") & (outcome_denominator=="per 100,000 pop at risk") & (scenario == scn) & year %in% as.character(1:30))
  df_tmp$outcome = gsub("_averted", "\naverted",df_tmp$outcome)
  df_tmp$outcome = gsub("_", " ",df_tmp$outcome)
  df_tmp$year=as.numeric((as.character(df_tmp$year)))
  
  p=ggplot(df_tmp, aes(x= year, y=value, ymin=CI_low, ymax=CI_high, fill=group, color=group, group=group)) +
    geom_line() +
    geom_ribbon(alpha=0.2, color=NA) +
    facet_grid(outcome~transmission_setting, scale="free_y") +
    xlab("time after introduction of CYD (years)") +
    ylab("per 100,000 pop at risk") +
    theme_bw() +
    scale_fill_manual(values=cbPalette) + scale_color_manual(values=cbPalette) 

  #print(p)
  ggsave(paste(path_figures,"VaccineEffectsInPopulation.",output,sep=""),p,dpi=300,units="cm",width=20, height=14)
  
  
  #change in serostatus in the targeted age group over time
  df_tmp=subset(df, outcome=="seropositive" & (scenario == scn))
  df_tmp$outcome = gsub("_averted", "\naverted",df_tmp$outcome)
  df_tmp$outcome = gsub("_", " ",df_tmp$outcome)
  df_tmp$year=as.numeric((as.character(df_tmp$year)))
  
  p=ggplot(df_tmp, aes(x= year, y=value, ymin=CI_low, ymax=CI_high, fill=group, color=group, group=group)) +
    geom_line() +
    geom_ribbon(alpha=0.2, color=NA) +
    facet_grid(.~transmission_setting, scale="free_y") +
    xlab("time after introduction of CYD (years)") +
    ylab("proportion seropositive\nat 9yrs of age") +
    theme_bw() +
    scale_fill_manual(values=cbPalette) + scale_color_manual(values=cbPalette) 
  
  ggsave(paste(path_figures,"VaccineEffectsOnSeropositivity.",output,sep=""),p,dpi=300,units="cm",width=20, height=8)
  
  

