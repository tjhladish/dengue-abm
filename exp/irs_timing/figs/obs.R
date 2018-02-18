require(data.table)

args <- commandArgs(trailingOnly = TRUE)

obs = fread(args[1],
  fill=T, col.names = sprintf("week_%02d",1:53)
)

obs.mn <- apply(obs[,1:52], 2, mean, na.rm=T) # ignore week 53 entry
# aside:
# obs.md <- apply(obs[,1:52], 2, median, na.rm=T)

wrapend = obs.mn[1]*(4.5/8) + obs.mn[52]*(3.5/8) # weight according to inverse distance from jan1/dec31
obs.mn = c(wrapend, obs.mn, wrapend)
norm.obs = obs.mn/max(obs.mn)

store <- data.table(
  doy=c(1, seq(3.5,365,7), 365),
  value = norm.obs,
  variable = "Cases (observed)",
  coverage = "reference",
  duration = "reference",
  durability = "reference",
  layer = "foreground"
)

saveRDS(store, args[2])
