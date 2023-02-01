library(DBI)
library(data.table)
library(lubridate)
library(ggplot2)
library(viridis)

if (interactive()) { setwd("~/documents/work/dengue/exp/tirs_study") }

.args <- if (interactive()) c(
  "sim_data_0.sqlite",
  "../../pop/merida-tirs/merida-tirs.sqlite"
) else commandArgs(trailingOnly = TRUE)

db_path <- .args[1]

serial <- sub("sim_data_([[:digit:]]*).sqlite", "\\1", db_path)
figPath <- paste0("./fig/sim_", serial, "_diagnostics")
dir.create(figPath, recursive = T)

# TIRS infections by location type
db <- dbConnect(RSQLite::SQLite(), db_path)
stmt <- paste0("attach database '", .args[2], "' as synthpop")
dbExecute(db, stmt)

# stmt <- "select ih.inf, ih.inf_owner_id, p.home_id as inf_owner_home, l.arm as inf_owner_arm, ih.inf_place_id, l.type as inf_place_type, ih.infected_time 
# from infection_history ih, synthpop.pop p, synthpop.loc l 
# where ih.inf_place_id = l.locid and ih.inf_owner_id = p.pid and p.home_id = l.locid and p.age >= 2 and p.age <= 15 and ih.infected_time >= 11099;"
# 
# inf_by_loc_history <- dbGetQuery(db, stmt)
# inf_by_loc_history <- setDT(inf_by_loc_history)

# stmt <- "select l.type, count(*) as num_inf 
# from infection_history ih, synthpop.loc l, synthpop.pop p 
# where ih.inf_place_id = l.locid and ih.inf_owner_id = p.pid and p.age >= 2 and p.age <= 15 and ih.infected_time >= 11099 
# group by l.type;"

stmt = "select ih.inf, ih.infected_time, ih.inf_place_id, l2.type, ih.inf_owner_id, p.home_id, l.arm, p.age 
from infection_history ih 
join synthpop.pop p on p.pid = ih.inf_owner_id 
join synthpop.loc l on l.locid  = p.home_id 
join synthpop.loc l2 on l2.locid = ih.inf_place_id 
where p.age >= 2 and p.age <= 15;" # and ih.infected_time >= 11099;"

inf_by_loc <- dbGetQuery(db, stmt)
inf_by_loc <- setDT(inf_by_loc)
inf_by_loc[arm == 1, arm_lab := "Control"]
inf_by_loc[arm == 2, arm_lab := "Treatment"]

dbDisconnect(db)

prop_TIRS_inf_by_loc <- ggplot(inf_by_loc[arm > 0, .(num_inf = .N), by = .(type, arm_lab)]) +
  geom_col(aes(x=factor(arm_lab), y=num_inf, fill=type), position = "fill") +
  labs(title = paste0("Serial ", serial), x = "TIRS trial arm", y = "% of infections in trial children", fill = "Location type") +
  theme(legend.position = "bottom")

ggsave(filename = paste0(figPath, "/TIRS_inf_by_loc.png"), plot = prop_TIRS_inf_by_loc, device = 'png', units = 'in', height = 6, width = 6, dpi = 300)