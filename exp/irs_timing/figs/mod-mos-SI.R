rm(list)
require(reshape2)
require(data.table)
require(lubridate)
require(ggplot2)

args <- c("~/Dropbox/who/fig1_data/mos.rds", "~/Dropbox/who/fig1_data/mos_mod.rds")
args <- commandArgs(trailingOnly = TRUE)

basemospop <- readRDS(args[1])
res <- readRDS(args[2])

doys <- 1:365-1 # need to make 0-364 to add a scale that starts w/ 1
dates <- lubridate::as_date(doys) # convert doys into dates for use w/ lubridate; default year 1970 is non-leap
month.len <- rle(month.abb[month(dates)])$lengths # get the month lengths
tick.locs <- c(0, cumsum(month.len)) + 1 # place ticks at 0, plus the cumsum of each month length

# want: labels falling in middle of month boundaries
name.locs <- head(tick.locs,-1) + month.len/2 #
mon.y <- -0.08

targetpoints <- merge(
  res[coverage==0.75 & efficacy == 0.8 & duration == 90],
  basemospop[layer=="foreground", value, by=doy],
on="doy")

targetpoints[, layer := "foreground"]
rainfall <- copy(targetpoints)[, layer := "background"][, multiplier := 1.0 ]
irs.dt <- targetpoints[,
  .(fillin=multiplier != 1, doy),
  by=.(start, durability)
]

pro.col <- "#0000ff"
rea.col <- "#664400"

p <- ggplot(
  rbind(targetpoints, rainfall)
) + facet_grid(durability ~ start, labeller = labeller(
  start=c(`148`="Proactive", `323`="Reactive"),
  durability=c(`30`="30 day durability",`90`="90 day durability",`150`="150 day durability")
)) +
  aes(x=doy, y=multiplier*value, linetype=layer, color="rainfall") +
  geom_line() +
  geom_rect(
    aes(xmin=doy, xmax=doy+1, ymin=0+2*mon.y, ymax=1-2*mon.y, y=NULL, linetype=NULL, color=NULL, fill=factor(start)),
    data=irs.dt[fillin==TRUE],
    alpha = 0.2
  ) +
  theme_minimal() + theme(
    panel.spacing=unit(1,"lines"),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(margin = margin(t=-12, unit="pt")),
    strip.text.y = element_text(angle=90)
  ) +
  coord_cartesian(xlim = c(1,365), ylim=c(0,1)) +
  scale_color_manual(values=c(rainfall="blue"), guide="none") +
  scale_fill_manual(values=c(`148`=pro.col,`323`=rea.col), guide="none") +
  scale_linetype_discrete(guide="none") +
  scale_y_continuous("Region-wide mosquito population",
                     breaks = (0:4)/4, minor_breaks = NULL,
                     labels = c(0,"","","","Max"),
                     expand=c(-mon.y*2,0)
  ) +
  scale_x_continuous("Time of year",
                     breaks = name.locs, minor_breaks = NULL,
                     labels = c("J","F","M","A","M","J","J","A","S","O","N","D"),
                     expand=c(0,0)
  ) +
  annotate("rect",
           xmin=head(tick.locs[c(TRUE, FALSE)],-1), # every other month start, except last
           xmax=tick.locs[c(FALSE, TRUE)], # every other month start, except first
           ymin=-Inf, ymax=Inf, # setting these to Inf ensures the rects extend full width + don't mess w/ scale
           alpha=0.05 # adjust for how dark; ALT: might also change color?
  )

#  +

ggsave(tail(args,1), width = unit(6.5,"in"), height = unit(5.5,"in"), dpi = 450, plot=p)
