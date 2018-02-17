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

foilvls <- c("70% FOI", "Baseline FOI", "130% FOI")
timelvls <- c("None", "Summer", "Winter")

ref.dt <- data.table(
  dbGetQuery(db, sprintf(
    'SELECT P.serial AS serial,
       CASE WHEN vector_control = 0 THEN "%s" WHEN timing = 152 THEN "%s" ELSE "%s" END AS intervention,
       CASE WHEN foi < 1.0          THEN "%s" WHEN foi = 1.0    THEN "%s" ELSE "%s" END AS foi
     FROM par P, job J
       WHERE P.serial = J.serial
         AND status = \'D\';',
      foilvls[1], foilvls[2], foilvls[3],
      timelvls[1], timelvls[2], timelvls[3]
  )), stringsAsFactors = T
)

ref.dt[,
  foi := factor(foi, levels = foilvls, ordered = T)
][,
  intervention := factor(intervention, levels = timelvls, ordered = T)  
]

## TODO check for .Rdata

daily.dt <- fread(
  log.path, nrows = 5, sep = " ", header = F,
  col.names = c("serial","day","incidence.intro","prevalence.intro","incidence.local","prevalence.local"),
  drop = "V1"
)

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
  facet_grid(foi ~ ., scales = "free") +
#  geom_ribbon(aes(ymin=min, ymax=max, alpha="full range", color=NULL)) +
  geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, alpha=views[1], color=NULL), data=irs.dt) +
  geom_ribbon(aes(x=day, ymin=lo, ymax=hi, alpha=views[2], color=NULL)) +
#  geom_line(aes(y=ave, linetype="average", alpha="line")) +
  geom_line(aes(x=day, y=md, linetype="median", alpha=views[3])) +
  scale_y_log10(name="# active infections") +
  scale_alpha_manual(name="measure", limits=views, values = c(0.1,0.2,1)) +
  guides(linetype="none")
