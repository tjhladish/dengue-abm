require(data.table)
require(ggplot2)
# inter-quartile plots

args <- c("effectiveness.rds","comboeff.rds")
args <- commandArgs(trailingOnly = TRUE)

effectiveness.dt <- readRDS(args[1])
comboeff.dt <- readRDS(args[2])

stats.dt <- effectiveness.dt[,{
  qs <- quantile(c.eff, probs = c(.25,.5,.75), na.rm = T)
  names(qs) <- c("IQlo","med","IQhi")
  as.list(qs)
}, keyby=.(vc, vac, vc_coverage, vac_mech, catchup, year)]

stats.dt[vc == 0, vc_coverage := 0]
stats.dt[, vac_mech := factor(ifelse(vac==1,c("cmdvi","trad")[vac_mech+1],"none"))]
stats.dt[, catchup := factor(c("none","catchup")[catchup+1])]

flabs <- labeller(
  catchup=c(none="None", catchup="Catchup"),
  vc_coverage=c(`0`="0%",`25`="25%",`50`="50%",`75`="75%")
)

gds <- function(
  order, title.position = "top",
  direction = "horizontal", ...
) guide_legend(
  title.position = title.position, direction = direction, order = order,
  label.position = "top", ...
)

p <- ggplot(stats.dt) +
  aes(x=year+1, y=med, ymin=IQlo, ymax=IQhi, fill=vac_mech) +
  facet_grid(catchup ~ vc_coverage, labeller = flabs) +
  geom_ribbon(alpha=.5) +
  geom_line(mapping=aes(color=vac_mech)) +
  theme_minimal() +
  scale_color_manual("Vaccine",
    labels=c(cmdvi="CMDVI",trad="Traditional", none="None"),
    values=c(cmdvi="blue", trad="green", none="black"),
    guide=gds(order=1)
  ) +
  scale_fill_manual("Vaccine",
    labels=c(cmdvi="CMDVI",trad="Traditional", none="None"),
    values=c(cmdvi="blue", trad="green", none="black"),
    guide=gds(order=1)
  ) +
  scale_x_continuous("Year", expand = c(0,0)) +
  scale_y_continuous("Cumulative Effectiveness", expand = c(0,0)) +
  coord_cartesian(ylim=c(-.2, 1)) +
  ggtitle("Vector Control Coverage") +
  theme(
    panel.spacing = unit(15,"pt"),
    plot.title = element_text(size=rel(.75), hjust = 0.5, vjust=0),
    legend.box = "horizontal",
    legend.position = c(0.5,0.5), legend.justification = c(0.5, 0.5),
    legend.text = element_text(size=rel(0.6)),
    legend.title = element_text(size=rel(0.7)),
    legend.title.align = 0.5
  )

ggsave(
  tail(args,1), p, device = "png",
  width = 7.5, height = 5, dpi = "retina", units = "in"
)