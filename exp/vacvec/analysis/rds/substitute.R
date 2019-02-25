suppressPackageStartupMessages({
	require(data.table)
})

args <- c("foi_intervention.rds","foi_test_intervention.rds")
args <- c("intervention.rds","test_intervention.rds")
args <- commandArgs(trailingOnly = TRUE)

orig <- readRDS(args[1])[vaccine != "cydtdv"]
updt <- readRDS(args[2])

saveRDS(setkeyv(rbind(orig,updt),key(orig)), tail(args, 1))
