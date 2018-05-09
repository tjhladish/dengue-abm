require(reshape2)
require(data.table)
require(lubridate)
require(ggplot2)

args <- c("~/Dropbox/who/fig1_data/mos.rds", "~/Dropbox/who/fig1_data/mod-mos.rds")
args <- commandArgs(trailingOnly = TRUE)

basemospop <- readRDS(args[1])
res <- readRDS(args[2])

ggsave(tail(args,1), width = unit(5.5,"in"), height = unit(7,"in"), dpi = 450, plot=
ggplot(
  merge(res[coverage==0.75 & efficacy == 0.8 & duration == 90], basemospop[layer=="foreground",value,by=doy], on="doy")
) + facet_grid(durability ~ start, labeller = labeller(
  start=c(`148`="Proactive", `323`="Reactive"),
  durability=c(`30`="30 day durability",`90`="90 day",`150`="150 day")
)) +
  aes(x=doy, y=multiplier*value
      #, alpha=factor(coverage), linetype=factor(duration), color=factor(efficacy)
  ) +
  geom_line() +
  annotate("line", x=basemospop[layer=="foreground",doy], y=basemospop[layer=="foreground", value], color="blue") +
  theme_minimal() +
  coord_cartesian(xlim = c(1,365), ylim=c(0,1)) +
  scale_y_continuous("Region-wide Mosquito Population") +
  scale_x_continuous("Day-of-Year")
)
