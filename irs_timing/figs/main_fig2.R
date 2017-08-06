rm(list=ls())

require(data.table)
require(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
# args <- c(paste0("~/Dropbox/who/fig1_data/stopping-",c("eff","sero"),".rds"),"~/Dropbox/who/fig1_data/fig2.png")

# should be 2 dimensions
# coverage (3 levels)
# stopping vs not (2 levels)
effectiveness.dt  <- readRDS(args[1])

# three cases
#  no intervention
#  75% coverage
#  75% coverage w/ stopping
seroprevalence.dt <- readRDS(args[2])

# two facets:
## effectiveness
## seroprevalence

eff.slice.dt <- effectiveness.dt[, .(value=q.med,measure="Annual effectiveness"), keyby=.(end_year, coverage, Year=year)]
sero.slice.dt <- seroprevalence.dt[, .(value=q.med,measure="Seroprevalence"), keyby=.(end_year, coverage, Year=year)]
sero.slice.dt[end_year == 0, end_year := 50L]

slice.dt <- rbind(eff.slice.dt, sero.slice.dt)

stopping_map <- c(`50`="Consistent annual IRS  ",`10`="IRS ends after 10 years")
slice.dt[,
  maintained := factor(
    stopping_map[as.character(end_year)],
    levels=stopping_map, ordered = T
  ) # use this factor to express different levels for linetype, establish labels
]

gg_combo_plot <- function(slice.dt, plot_years=20) {
  dt <- slice.dt[Year <= plot_years]
  space.dt <- data.table(
    end_year=rep(50,4),
    coverage=rep(0,4),
    Year=rep(0,4),
    value = c(1,-4,1,.4), # boundaries
    measure = c("Annual effectiveness","Annual effectiveness","Seroprevalence","Seroprevalence"),
    maintained=rep(stopping_map[1],4)
  )
  labels.dt <- data.table(
    Year = c(1,1), value = c(-3,.90), measure = c("Annual effectiveness","Seroprevalence"),
    lab = c("a","b"), maintained = rep(stopping_map[1],2)
  )
  
  ggplot(dt,
         aes(
           x = Year, y = value,
           color = factor(coverage, levels=c(75,50,25,0), ordered = T),
           # factor to get discrete scale
           # reverse levels to get black -> grey :: 75 -> 25
           linetype = factor(maintained)
         )
  ) + facet_grid(measure ~ ., scales = "free_y", switch = "y") + 
    geom_line(size=1) + geom_blank(data=space.dt) +
    geom_text(aes(label=lab, color=NULL), labels.dt, size=20, show.legend = F) +
    theme_minimal() + theme(
      text = element_text(size = 25),
      panel.background = element_rect(fill="grey99", color=NA),
      
      legend.box.background = element_rect(fill=alpha('white',.6), color=NA),
      legend.margin = margin(),
      legend.text = element_text(size=rel(0.8)),
      legend.key.width = unit(1, "cm"),
      legend.key.height = unit(0.9, "cm"),
      legend.title = element_text(face="bold", size=rel(0.8)),
      legend.direction = "horizontal",
      legend.position = c(0.5, 0.45),
      
      axis.title.x = element_text(size=rel(1)),
      
      strip.placement = "outside",
      strip.text = element_text(size=rel(1)),
      panel.grid.minor = element_blank(),
      panel.spacing.y = unit(1,"cm")
      
    ) +
    scale_y_continuous(name="") + scale_x_continuous(name="Years since introducing IRS") +
    scale_color_grey(name="IRS Coverage", labels=function(l) { 
      c(sprintf("%s%%  ",head(l,-1)),"0% (baseline)")
    }) +
    scale_linetype(name="Intervention") +
    guides(
      linetype=guide_legend(title.position = "top", order = 1),
      color=guide_legend(title.position = "top", order = 2)
    )
}

res = 300
mag = 0.85*res/72
png(args[3], width = 1000*mag, height = 1000*mag, units = "px", res=res)
gg_combo_plot(slice.dt)
dev.off()

# plot_effectiveness_over_time = function(eff.dt, plot_years=20) {
#   slice.dt <- eff.dt[year < plot_years, q.med, keyby=.(end_year, coverage, year)]
#   # keyby sets order ascending by end_year, coverage, year
#   # increasing 'end_year' -> all the 10s before the 50s
#   # next, increasing 'coverage' -> within 10/50, 25 -> 50 -> 75
#   # last, increasing year w/in each scenario
#   # this is going to lead to drawing 10 first, then 50
# 
#   # build matrix column-by-column from slice.dt
#   ys <- matrix(slice.dt$q.med, byrow=F, nrow=plot_years)
#   reds=c('lightgrey','grey','black') # increasing tone w/ increasing coverage
#   .lwd = 1 # only one line weight
#   lty10=3; lty50=1 # ltys by 10 vs 50
#   ylim_ = c(-4.4, 1.4)
#   matplot(y=ys, type='l',
#     lwd=.lwd, lty=c(rep(lty10,3), rep(lty50,3)),
#     col=c(reds, reds),
#     ylab='Overall effectiveness', xlab='', axes=F, ylim = ylim_
#   )
# 
#   ## tjh: TODO, if you want legend-y stuff sorted
#   text(1, ylim_[2]-(0.025*(ylim_[2]-ylim_[1])), cex=2, font=2,labels='a')
#   box()
#   axis(1, labels=F)
#   axis(1, lwd = 0, line = -0.3)
#   axis(2)
# 
#   abline(h=0,lty=2)
#   legend('bottomleft', legend=c('75% coverage','50% coverage','25% coverage','75%, ended year 10','50%, ended year 10','25%, ended year 10'), lwd=rev(.lwd), lty=1:3, bty='n', col=c(rev(reds),rep('#999999',3)),cex=0.8)
#   return(medians)
# }
