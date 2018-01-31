
require(data.table)
require(ggplot2)
require(lubridate)

args <- commandArgs(trailingOnly = TRUE)
# args <- paste0("~/Dropbox/who/fig1_data/",c(paste0(c("eip","R0","mos"),".rds"), "coverage_illustration.png"))

eip.dt <- readRDS(args[1])
R0.dt <- readRDS(args[2])
mos.dt <- readRDS(args[3])

y.lim <- 100 # alt: ylim <- 1

optimal_start_date <- 148
covered_households <- 0.75*y.lim

doys <- 1:365-1 # need to make 0-364 to add a scale that starts w/ 1
dates <- lubridate::as_date(doys) # convert doys into dates for use w/ lubridate; default year 1970 is non-leap
month.len <- rle(month.abb[month(dates)])$lengths # get the month lengths
tick.locs <- c(0, cumsum(month.len)) + 1 # place ticks at 0, plus the cumsum of each month length

# want: labels falling in middle of month boundaries
name.locs <- head(tick.locs,-1) + month.len/2 #

intervention_examples.dt <- data.table(
  coverage = c(
    rep((covered_households/365)*90, 2),
    c(0, 0, rep(covered_households, 2), 0, 0),
    c(0, 0, covered_households, 0, 0)
  ),
  doy=c(
    c(1, 365),
    c(1, rep(optimal_start_date+45,2), rep(optimal_start_date+45+90,2), 365),
    c(1, optimal_start_date, optimal_start_date+90, optimal_start_date+180, 365)
  ),
  intervention=factor(
    c(rep("continuous", 2), rep("instant", 6), rep("90", 5)),
    levels=c("instant","90","continuous"),
    ordered = T
  )
)

y.exp <- 0.07
mon.y <- -y.exp/2*y.lim

## TODO: add seasonal factors in background?

ggsave(tail(args,1), plot=ggplot(intervention_examples.dt) + theme_minimal() + theme(
  # keep on y.major grid lines
  panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
  # legend size / locating:
  legend.key.width = unit(45/365, "npc"),
  # npc units are 0-1 on the plot space; width is x dim
  #   which is 365 long, so this means make the key 45 days wide
  legend.justification = c("left", "top"), # position which corner of the legend; includes margins
  legend.margin = margin(0,0,0,0), # zero out legend margins, so position off minimum bound of key/labels
  legend.position = c(0, 0.95), # npc coordinates, so roughly put corner == day 1, coverage = 0.95*y.lim
  legend.title = element_blank() # don't title the legend
) +
  aes(x=doy, y=coverage, linetype=intervention) +
  geom_line() +
  # geom_line() + # TODO add seasonal curves to background?
  coord_cartesian(ylim=c(0,y.lim)) +
  scale_y_continuous("Coverage (%)", expand = c(y.exp, 0)) + # y.exp provides space to put month labels
  scale_x_continuous("Julian Day", breaks = NULL, expand=c(0,0)) +
  scale_linetype_manual(
    values=c(instant=2, continuous=3, `90`=1),
    labels=c(instant='all covered houses treated in 1 day', continuous='covered houses treated across 365 days', `90`='covered houses treated across 90 days')
  ) +
  annotate("text",
    x = name.locs, y = mon.y,
    label = c("J","F","M","A","M","J","J","A","S","O","N","D"),
    size = rel(4)
  ) + annotate("rect",
    xmin=head(tick.locs[c(TRUE, FALSE)],-1), # every other month start, except last
    xmax=tick.locs[c(FALSE, TRUE)], # every other month start, except first
    ymin=-Inf, ymax=Inf, # setting these to Inf ensures the rects extend full width + don't mess w/ scale
    alpha=0.05 # adjust for how dark; ALT: might also change color?
  ),
  width=unit(6,"in"), height=unit(4,"in"), dpi=300
)
