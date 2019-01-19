require(data.table)
require(ggplot2)
# inter-quartile plots

source("projref.R")

args <- c("effectiveness.rds")
args <- commandArgs(trailingOnly = TRUE)

effectiveness.dt <- readRDS(args[1])

stats.dt <- effectiveness.dt[,{
  qs <- quantile(eff, probs = c(.25,.5,.75), na.rm = T)
  names(qs) <- c("IQlo","med","IQhi")
  as.list(qs)
}, keyby=.(vc, vac, vc_coverage, vaccine, catchup, year)]

flabs <- labeller(
  catchup=c(none="None", catchup="Catchup"),
  vc_coverage=c(`0`="0%",`25`="25%",`50`="50%",`75`="75%")
)

p <- ggplot(stats.dt) +
  aes(x=year+1, y=med, ymin=IQlo, ymax=IQhi, fill=vaccine) +
  facet_grid(catchup ~ vc_coverage, labeller = flabs) +
  geom_ribbon(alpha=.5) + geom_line(mapping=aes(color=vaccine), show.legend = F) +
  theme_minimal() +
  scale_x_continuous("Year", expand = c(0,0)) +
  scale_y_continuous("Effectiveness") +
  scale_color_manual("Vaccine", values=vac_cols) +
  scale_fill_manual("Vaccine", values=vac_cols) +
  ggtitle("Vector Control Coverage") +
  theme(
    panel.spacing = unit(15,"pt"),
    plot.title = element_text(size=rel(.75), hjust = 0.5, vjust=0)
  )

ggsave(
  tail(args,1), p, device = "png",
  width = 7.5, height = 5, dpi = "retina", units = "in"
)