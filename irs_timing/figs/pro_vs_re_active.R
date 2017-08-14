rm(list=ls())

args <- commandArgs(trailingOnly = T)
# args <- c("~/Dropbox/who/fig1_data/campaign-timing-%.rds", "~/Dropbox/who/fig1_data/pro_vs_re_active.png")

require(data.table)
require(ggplot2)

plot.dt <- readRDS(args[1])
timelvls <- c("None", "Proactive", "Reactive")

start.day <- 125; end.day <- 900
ymin <- 1; ymax <- 1e4; reps <- 5
starts <- c(seq(from=152,by=365,length.out=reps),seq(from=305,by=365,length.out=reps))
irs.dt <- data.table(
  xmin=starts,
  xmax=starts+89+90,
  ymin=ymin, ymax=ymax,
  intervention = factor(c(rep(timelvls[2],reps),rep(timelvls[3],reps)), levels = timelvls, ordered = T)
)
views = c("IRS window","Half range","Median")
view <- function(i) {
  views = c("IRS window","Half range","Median")
  factor(views[i], levels = views, ordered = T)
}

# maybe use this:
# https://stackoverflow.com/questions/36941197/overall-label-for-facets?noredirect=1&lq=1

# plotter <- function(view, ymin=ymin, ymax=ymax) ggplot(plot.dt[variable == view],
#   aes(color=intervention,fill=intervention)
# ) + theme_minimal() +
#   coord_cartesian(xlim=c(start.day,end.day), ylim=c(ymin,ymax)) +
#   facet_grid(foi ~ ., scales = "free") +
#   #  geom_ribbon(aes(ymin=min, ymax=max, alpha="full range", color=NULL)) +
#   geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, alpha=views[1], color=NULL), data=irs.dt) +
#   geom_ribbon(aes(x=day, ymin=lo, ymax=hi, alpha=views[2], color=NULL)) +
#   #  geom_line(aes(y=ave, linetype="average", alpha="line")) +
#   geom_line(aes(x=day, y=md, linetype="median", alpha=views[3])) +
# #  scale_y_log10(name="# Active Infections (log scale)") +
#   scale_alpha_manual(name="measure", limits=views, values = c(0.1,0.2,1)) +
#   guides(linetype="none")

pro.col <- "#0000ff"
rea.col <- "#664400"

png(args[2],height=1200,width=1000,res=180)
ggplot(plot.dt,
       aes(fill=intervention)
) + theme_minimal() + theme(strip.text.y = element_text(angle=90)) +
  coord_cartesian(xlim=c(125,850), ylim=c(ymin,ymax)) +
  facet_grid(foi ~ ., scales = "free") +
  #  geom_ribbon(aes(ymin=min, ymax=max, alpha="full range", color=NULL)) +
  geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, alpha=views[1]), data=irs.dt) +
  geom_text(aes(x=(xmax+xmin)/2, y=ifelse(intervention == "Proactive",ymax*.6,ymin/.6), label=intervention, color=intervention), data=irs.dt) +
  geom_ribbon(aes(x=day, ymin=lo, ymax=hi, alpha=views[2])) +
  #  geom_line(aes(y=ave, linetype="average", alpha="line")) +
  geom_line(aes(x=day, y=md, linetype="median", alpha=views[3], color=intervention)) +
  scale_y_log10(name="# Active infections") +
  scale_alpha_manual(name="Measure", limits=views, values = c(0.1,0.2,1)) +
  scale_color_manual(name="Intervention", values=c(None='black',Proactive=pro.col,Reactive=rea.col)) +
  scale_fill_manual(name="Intervention", values=c(None='black',Proactive=pro.col,Reactive=rea.col)) +
  scale_x_continuous(name="Days since start of first intervention calendar year") +
  guides(linetype="none")
dev.off()