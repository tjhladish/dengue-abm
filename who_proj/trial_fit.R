## read in sim outputs
## out format: Groups Trial    Description    Arm       Strata        mean       lower       upper

require(RSQLite)

getMetrics <- function(sqlite) {
  db <- RSQLite::dbConnect(RSQLite::SQLite(), sqlite)
  result <- data.table(RSQLite::dbGetQuery(db, "select seed, * from metrics m, parameters p where m.serial = p.serial;"))
  RSQLite::dbDisconnect(db)
  result
}

init <- getMetrics("~/Downloads/cyd14_match.sqlite")[beta_multiplier == 4] #, select=-c(serial, seed, seed, mild_EF, severe_EF, base_path, sec_sev, pss_ratio, exp_coef, num_mos, beta, beta_multiplier)

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
  se = sd(value)/sqrt(.N)
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

saveRDS(subset(wider,select=-variable)[,Groups:="longini"], "~/Dropbox/CMDVI/Phase II analysis/Data/UF-Longini/longini-fit.rds")
