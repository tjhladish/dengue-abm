suppressPackageStartupMessages({
  require(data.table)
})

warnnonunique <- function(var, variable, collapse = median) {
  if (length(unique(var)) != 1) warning(sprintf("non unique %s", unique(variable)))
  collapse(var)
}

# debugging args for interactive use
args <- c("projref.rda", "effstats.rds", "base_single_intervention.rds")

args <- commandArgs(trailingOnly = TRUE)

# load the reference digests
load(args[1])
effstats.dt <- readRDS(args[2])
tar <- tail(args, 1)

vac.eff <- effstats.dt[variable == "vac.eff", .(
    value = warnnonunique(med, variable),
    vc_coverage = 0,
    scenario = trans_scnario(0, 1)
  ), keyby=.(
    vaccine, catchup = ifelse(catchup=="vc+vac","vac-only","routine"), year
    #, measure = trans_meas(gsub("vac\\.","",variable))
  )
]

vec.eff <- cbind(effstats.dt[variable == "vec.eff", .(
    value = warnnonunique(med, variable),
    scenario = trans_scnario(1, 0)
  ), keyby=.(
    vc_coverage, year
    #, measure = trans_meas(gsub("vec\\.","",variable))
  )
], reference.scenario[,.(vaccine, catchup)])

plot.dt <- rbind(vec.eff, vac.eff)

saveRDS(plot.dt, tar)