suppressPackageStartupMessages({
	require(data.table)
})

args <- c("foi_intervention.rds","foi_test_intervention.rds")
args <- commandArgs(trailingOnly = TRUE)

orig <- readRDS(args[1])[scenario == "vc"]
updt <- readRDS(args[2])

saveRDS(setkeyv(rbind(orig,updt),key(orig)), tail(args, 1))