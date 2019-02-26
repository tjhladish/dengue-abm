suppressPackageStartupMessages({
	require(data.table)
})

args <- c("foi_intervention.rds","foi_test_intervention.rds")
args <- c("intervention.rds","test_intervention.rds")
args <- commandArgs(trailingOnly = TRUE)

orig <- readRDS(args[1])
keep <- orig[vaccine != "cydtdv"]
leftover <- orig[vaccine == "cydtdv"]
updt <- readRDS(args[2])

saveRDS(setkeyv(rbind(keep, updt),key(orig)), tail(args, 1))
saveRDS(leftover, gsub("mrg","sub",tail(args, 1)))