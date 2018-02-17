#!/apps/R/3.2.0/bin/Rscript

require(data.table)
path_to_event_logs <- commandArgs(trailingOnly = TRUE)[1]

filenames <- dir(path_to_event_logs, pattern = "daily")

thing <- sapply(filenames, function(filename, p_severe_nonhosp, p_death_nonhosp, p_severe_hosp, p_death_hosp) { 
  whorun.dt <- fread(filename)[,
    prior := nchar(gsub("-","", hist))-1
  ][,
    infection := factor(
      ifelse(prior == 0, "primary", ifelse(prior == 1, "secondary", "post-secondary")),
      levels = c("primary","secondary","post-secondary")
    )
  ][,
    year := floor(day / 365)
  ][,
    outcome_base := factor(
      ifelse(!case, "asymptomatic", ifelse(!hosp, "clinical.case", "hospitalised")),
      levels = c("asymptomatic", "clinical.case", "hospitalised")
    )
  ][,
    outcome_draw := runif(.N)
  ]
  
  whorun.dt[ outcome_base == "asymptomatic",  outcome_adv := factor("none", levels=c("none", "severe", "death"))]
  whorun.dt[ outcome_base == "clinical.case", outcome_adv := factor(ifelse(outcome_draw < p_severe_nonhosp, "none", ifelse(outcome_draw < p_death_nonhosp, "severe", "death")), levels=c("none", "severe", "death"))]
  whorun.dt[ outcome_base == "hospitalised",  outcome_adv := factor(ifelse(outcome_draw < p_severe_hosp,    "none", ifelse(outcome_draw < p_death_hosp,    "severe", "death")), levels=c("none", "severe", "death"))]
  
  whorun.dt[,
    .N, keyby=list(year, serotype=sero+1, age, infection, outcome_base, outcome_adv)
  ]
  save(whorun.dt, file=sub("daily", "digest", filename))
}, p_severe_nonhosp = (1-1/54), p_death_nonhosp = (1-1/54*0.0103), p_severe_hosp = 1/2, p_death_hosp = (1-1/2*0.0103), simplify = F, USE.NAMES = F)

# sqlite_file <- "who-toy.sqlite"
# otherthing <- abc.import(sqlite_file, list(parameters=NULL))
# filenames <- paste0(path_to_event_logs,"daily.1870348530.-1379691560")

# slice <- melt(
#   otherthing$metrics[otherthing$parameters][vac==0,.SD,.SDcols=c("serial", grep("tc", names(otherthing$metrics), value = T))], id.var ="serial"
# )[,
#   serotype := sub("s(\\d).+","\\1", variable)
# ][,
#   year := as.integer(sub("s\\dtc(\\d+)","\\1", variable))
# ]
# 
# qplot(year, cases, 
#   data = slice[,list(value=sum(value)),keyby=list(year, serial)][,list(cases = mean(value)/18.2),keyby=year], geom="line", ylim = c(0,NA)) +
#   theme_bw() + scale_x_discrete() + ylab("cases per 100k")
# 
# whorun.dt <- fread("../../../../dengue_tmp_data/who_toy_data/scratch/lfs/thladish/who_output/daily.60086540.-1563855759")
# 
# whorun.dt[, prior := nchar(gsub("-","", hist))-1 ]
# whorun.dt[, infection := ifelse(prior == 0, "primary", ifelse(prior == 1, "secondary", "post-secondary")) ]
# whorun.dt[, year := floor(day / 365) ]
# 
# thing <- whorun.dt[,.N,keyby=list(age, sero, infection, year)]
# 
# require(ggplot2)
# 
# ggplot(thing) + theme_bw() +
#   facet_grid(year ~ .) + aes(x=age, y=N, fill=factor(sero)) + geom_bar(stat="identity")
