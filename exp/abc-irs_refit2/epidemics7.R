suppressPackageStartupMessages({
  require(data.table)
  require(ggplot2)
  require(cowplot)
})

.args <- c("~/Downloads/posterior_annual_incidence.out", "~/Dropbox/who/epidemics7_in.txt", "~/Dropbox/who/epidemics7.png")
.args <- commandArgs(trailingOnly = T)

posterior.dt <- fread(.args[1])
pop100ks <- 18.2
ref <- posterior.dt[,{
  qs <- quantile(case*mild_rf + dss_rf*dss, probs = c(0.025,0.5,0.975))/pop100ks
  .(lo=qs[1], md=qs[2], hi=qs[3])
},by=year+1879]

obs.dt <- fread(.args[2], col.names = c("year","cases"))

periods <- data.table(
  start = c(1879, 1956, 1979, 2015),
  end = c(1956, 1979, 2015, 2030),
  label = c("Priming", "DDT", "Fitting", "Forecast")
)
periods[, mid := (start+end)/2 ]

p <- ggplot(ref) + aes(year) +
  geom_vline(aes(xintercept=brks), data.table(brks=c(1955, 1979, 2015)), color="grey") +
  geom_ribbon(aes(ymin=lo, ymax=hi, fill="sim"), alpha=0.3) +
  geom_line(aes(y=md, color="sim")) +
  geom_line(aes(y=cases, color="obs"), obs.dt) +
  geom_text(aes(y=1150, x=mid, label=label), periods, fontface="bold", size=4.5) +
  scale_color_manual("",
    labels = c(sim="Simulated Median (95% IR)", obs="Observed"),
    values = c(sim="firebrick", obs="black"),
    aesthetics = c("color","fill"), guide = guide_legend(
      override.aes = list(fill=c(obs=NA, sim="firebrick"))
    )
  ) +
  scale_y_continuous("Annual Reported Cases Per 100k People", expand = c(0,0)) +
  scale_x_continuous(
    "Year",
    expand = c(0,0),
    breaks = seq(1880, 2020, by=20)
  ) +
  coord_cartesian(ylim = c(0,1200)) +
  theme_minimal() + theme(
    title = element_text(size=14),
    legend.text = element_text(size=rel(1.25)),
    axis.title = element_text(size=rel(1)),
    axis.text = element_text(size=rel(1.25)),
    axis.text.x = element_text(vjust=0),
    legend.position = c(0.25, 0.8),
    panel.grid = element_blank()
  )

save_plot(tail(.args, 1), p, base_width = 11, base_height = 5)
