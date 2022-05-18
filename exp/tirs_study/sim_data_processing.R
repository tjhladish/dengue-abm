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

stmt <- "select l.type, count(*) as num_inf 
from infection_history ih, synthpop.loc l, synthpop.pop p 
where ih.inf_place_id = l.locid and ih.inf_owner_id = p.pid and p.age >= 2 and p.age <= 15 and ih.infected_time >= 11099 
group by l.type;"

inf_by_loc_ct <- dbGetQuery(db, stmt)
inf_by_loc_ct <- setDT(inf_by_loc_ct)

dbDisconnect(db)

prop_TIRS_inf_by_loc <- ggplot(inf_by_loc_ct) +
  geom_col(aes(x=paste0("Simulation serial ", serial), y=num_inf, fill=type), position = "fill") +
  labs(x = "", y = "% of infections in trial children")

ggsave(filename = paste0(figPath, "/TIRS_inf_by_loc.png"), plot = prop_TIRS_inf_by_loc, device = 'png', units = 'in', height = 12, width = 6, dpi = 300)