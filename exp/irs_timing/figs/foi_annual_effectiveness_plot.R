rm(list=ls())

require(data.table)
require(ggplot2)
require(gridExtra)
require(grid)

args <- commandArgs(trailingOnly = T)
# args <- c("~/Dropbox/who/fig1_data/foi-baseline.cases.rds", "~/Dropbox/who/fig1_data/foi-interventions.cases.rds", "~/Dropbox/who/fig1_data/foi-averted.rds", "~/Dropbox/who/fig1_data/foi-eff.rds", "~/Dropbox/who/fig1_data/4panel_foi_20yr.png")

plot_years <- as.integer(gsub(".+_(\\d+)yr.png","\\1", args[5]))

base.dt <- readRDS(args[1])
bkeys <- grep("(case|particle)", names(base.dt), invert=T, value=T)

baseline.cases.dt <- base.dt[year < plot_years,.(
  cases.md = median(cases), cum.cases.md = median(cum.cases)
), by=bkeys ]

i.dt <- readRDS(args[2])
ikeys <- grep("(case|particle)", names(i.dt), invert=T, value=T)
plotkeys <- grep("year",ikeys,invert = T, value = T)

interventions.cases.dt <- i.dt[year < plot_years,.(
  cases.md = median(cases), cum.cases.md = median(cum.cases)
), by=ikeys]

averted.cases.dt <- readRDS(args[3])[year < plot_years]
eff.dt <- readRDS(args[4])[year < plot_years]

pop_size = 18.2 # in 100 thousands

################
# ggplot version

baseline.cases.dt[,
  timing := factor("baseline", levels=c("proactive","reactive","baseline"), ordered=T)
]

factorize <- function(dt) dt[,
  timing := factor(
    ifelse(timing == 152,"proactive","reactive"),
    levels = c("proactive","reactive","baseline"),
    ordered = T
  )
]

factorize(interventions.cases.dt)
factorize(eff.dt)
factorize(averted.cases.dt)

plotdata.dt <- rbind(baseline.cases.dt, interventions.cases.dt)

basep <- ggplot(averted.cases.dt, aes(
  x=year,
  color=factor(foi), size=factor(foi),
  linetype=timing, group=interaction(foi, timing)
)) + scale_color_manual(
  "Mosquito\nPopulation",
  labels = c(
    `0.7`=expression("70% M"[peak]),
    `1`="baseline",
    `1.3`=expression("130% M"[peak])
  ),
  values = c(`0.7`="blue",`1`="black",`1.3`="red")
) + scale_size_manual(
  "Mosquito\nPopulation",
  labels = c(
    `0.7`=expression("70% M"[peak]),
    `1`="baseline",
    `1.3`=expression("130% M"[peak])
  ),
  values = c(`0.7`=1,`1`=10/7,`1.3`=13/7)
) + scale_linetype_discrete(
  "Intervention\nTiming",
  labels=function(ls) gsub("^(.)", "\\U\\1", ls, perl=T)
) + scale_x_continuous("Year") + theme_classic() + theme(
  legend.box = "horizontal"
)

casep <- basep + geom_line(aes(y=cases.md/pop_size), data=plotdata.dt) + 
  scale_y_continuous("Cases per 100,000 people")
avertedp <- basep+geom_line(aes(y=cum.averted.md/(pop_size*1000))) +
  scale_y_continuous("Cumulative cases averted per 100 people")
effp <- basep+geom_line(aes(y=q.med), data=eff.dt) +
  scale_y_continuous("Annual Effectiveness") + coord_cartesian(ylim=c(0,1))
ceffp <- basep+geom_line(aes(y=cum.eff.md)) +
  scale_y_continuous("Cumulative Effectiveness") + coord_cartesian(ylim=c(0,1))

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}

png(args[5], width=2000, height=1720, res=225)
grid_arrange_shared_legend(
  casep, avertedp, effp, ceffp, ncol=2, nrow=2
)
dev.off()
