
args <- commandArgs(trailingOnly = TRUE)
# args <- paste0("~/Dropbox/who/fig1_data/stopping-",c("baseline","interventions","sero-SI"),".rds")

require(data.table)

baseline.dt <- readRDS(args[1])
interventions.dt <- readRDS(args[2])

extract <- function(x, pre) {
  qs <- quantile(x, probs = (0:4)/4)
  res <- list(qs[1], qs[2], qs[3], qs[4], qs[5], mean(x))
  names(res) <- paste(pre,c("min","lo","med","hi","max","mn"),sep=".")
  res
}

block <- expression({
  c(extract(imm0,"imm0"), extract(imm1,"imm1"), extract(imm2,"imm2"), extract(imm3,"imm3"), extract(imm4,"imm4"))
})

sero.dt <- rbind(
  interventions.dt[,
    eval(block),
    by=.(year, coverage, duration, end_year)
  ],
  baseline.dt[,
    eval(block),
    by=.(year)
  ][, c("coverage", "duration", "end_year") := .(0,0,0)]
)

saveRDS(sero.dt, args[3])
