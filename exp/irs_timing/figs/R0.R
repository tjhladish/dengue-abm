
require(reshape2)
require(data.table)
require(RSQLite)

args <- commandArgs(trailingOnly = TRUE)

drv = dbDriver("SQLite")
db  = dbConnect(drv, args[1], flags=SQLITE_RO)

all_r0 = dbGetQuery(db, 'select realization, m.* from met m, par p where m.serial=p.serial')
stat_r0 <- t(apply(all_r0[,-c(1,2)],2,function(col){
  # ignoring first two columns (realization / serial), get the quantiles of each col (day-of-year R0 observation)
  # then transpose, so that each row == doy observation
  c(quantile(col,probs = c(0,0.025,.25,.5,.75,0.975,1)), mean=mean(col))
}))

dbDisconnect(db)

# add column for smoothed mean
# replicate data 3x, and take center third
stat_r0 <- cbind(stat_r0, smooth=smooth.spline(rep(stat_r0[,"mean"],3),spar=0.5)$y[1:dim(stat_r0)[1]+dim(stat_r0)[1]])

store <- melt(
  data.table(
    doy = as.integer(gsub("r","",rownames(stat_r0)))+1,
    q = stat_r0
  ), id.vars = "doy", variable.name = "layer"
)[layer %in% c("q.mean", "q.smooth")]

norm <- store[layer == "q.smooth", max(value)]
store <- rbind(store[,
  .(doy, value=value/norm, variable="R0",
    coverage="reference", duration="reference", durability="reference",
    layer = ifelse(layer == "q.smooth","foreground","background")
)],
  data.table(
    doy=1, value=1/norm, variable="R0",
    coverage="reference", duration="365", durability="reference",
    layer="foreground"
  )
)

saveRDS(store, args[2])

# aside:
# ggplot(mapping = aes(x=doy)) + theme_minimal() +
#  geom_line(aes(y=value, linetype="smooth"), R0.dt[variable == "q.smooth"]) +
#  geom_line(aes(y=value, linetype="mean"), R0.dt[variable == "q.mean"]) +
#  geom_step(aes(y=value, linetype="median"), R0.dt[variable == "q.50%"]) +
#  geom_ribbon(
#    aes(ymin = `q.25%`, ymax=`q.75%`),
#    recast(R0.dt[variable %in% c("q.25%","q.75%")], doy ~ variable, measure.var = "value"),
#    alpha = 0.25
#  ) +
#  geom_ribbon(
#    aes(ymin = `q.97.5%`, ymax=`q.2.5%`),
#    recast(R0.dt[variable %in% c("q.2.5%","q.97.5%")], doy ~ variable, measure.var = "value"),
#    alpha = 0.1
#  ) +
#  geom_ribbon(
#    aes(ymin = `q.0%`, ymax=`q.100%`),
#    recast(R0.dt[variable %in% c("q.100%","q.0%")], doy ~ variable, measure.var = "value"),
#    alpha = 0.05
#  ) +
#  coord_cartesian(ylim = c(1,100)) + scale_y_log10(name="R_0")
