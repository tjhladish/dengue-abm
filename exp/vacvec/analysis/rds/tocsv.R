suppressPackageStartupMessages({
  require(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

srcrds <- readRDS(args[1])

fwrite(srcrds, tail(args,1), row.names = F)
