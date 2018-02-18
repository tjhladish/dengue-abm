
require(reshape2)
require(data.table)

args <- commandArgs(trailingOnly = TRUE)

raw <- fread(args[1], col.names = c("EIP","doy"))
reshaped <- melt(raw, id.vars = "doy")
mod <- reshaped[,.(
  doy, value=value/max(value), variable,
  coverage = "reference", duration = "reference", durability = "reference", layer = "foreground"
)]

saveRDS(mod, args[2])
