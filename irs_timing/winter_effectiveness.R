rm(list=ls())
require("RSQLite")
require("data.table")
drv = dbDriver("SQLite")
path <- sprintf("%s%s",.dbox.path,"who")
db.path <- sprintf("%s/%s",path,"irs_timing-summer_winter-foi.sqlite")
log.path <- sprintf("%s/%s",path,"daily-irs_refit-simple-nofilenames-uniq.out")

db = dbConnect(drv, db.path)

# campaign_duration is always 1
# coverage is always 75%

ref.dt <- data.table(
  dbGetQuery(db,
    'SELECT P.serial AS serial,
       vector_control,
       CASE timing WHEN 152 THEN "summer" ELSE "winter" END AS intervention,
       CASE WHEN foi=1.0 THEN "fit" WHEN foi<1.0 THEN "low" ELSE "high" END AS foi
     FROM par P, job J
       WHERE P.serial = J.serial
         AND status = \'D\';'
  )
)

ref.dt[, foi := factor(foi, levels = c("low","fit","high"), ordered = T)]
ref.dt[vector_control == 0, intervention := "none" ]
ref.dt$vector_control <- NULL

daily.dt <- fread(log.path,sep = " ", header = F)
setnames(daily.dt, c("V1","V2","V3","V4","V5","V6","V7"),
         c("V1","serial","day","incidence.intro","prevalence.intro","incidence.local","prevalence.local"))
daily.dt$V1 <- NULL

require(reshape2)

plot.dt <- setkey(melt.data.table(daily.dt, id.vars=c("serial","day")),serial,day)[ref.dt, on="serial"][,{
  mn <- mean(value)
  ps <- quantile(value, probs = c(0,.25,.5,.75, 1))
  list(ave=mn, min=ps[1], lo=ps[2], md=ps[3], hi=ps[4], max=ps[5])
}, keyby=list(intervention, foi, variable, day)
]
rm(daily.dt)
plot.dt[, day    := day - min(day) ]


require(ggplot2)

ymin <- 1; ymax <- 1e4; reps <- 5
starts <- c(seq(from=152,by=365,length.out=reps),seq(from=305,by=365,length.out=reps))
irs.dt <- data.table(
  xmin=starts,
  xmax=starts+89,
  ymin=ymin, ymax=ymax,
  intervention = c(rep("summer",reps),rep("winter",reps))
)
views = c("IRS window","half range","median")
view <- function(i) {
  views = c("IRS window","half range","median")
  factor(views[i], levels = views, ordered = T)
}

ggplot(plot.dt[variable == "prevalence.local"],
  aes(color=intervention,fill=intervention)
) + theme_minimal() +
  coord_cartesian(xlim=c(125,775), ylim=c(ymin,ymax)) +
  facet_grid(. ~ foi, scales = "free") +
#  geom_ribbon(aes(ymin=min, ymax=max, alpha="full range", color=NULL)) +
  geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, alpha=views[1], color=NULL), data=irs.dt) +
  geom_ribbon(aes(x=day, ymin=lo, ymax=hi, alpha=views[2], color=NULL)) +
#  geom_line(aes(y=ave, linetype="average", alpha="line")) +
  geom_line(aes(x=day, y=md, linetype="median", alpha=views[3])) +
  scale_y_log10(name="# active infections") +
  scale_alpha_manual(limits=views, values = c(0.1,0.2,1))
