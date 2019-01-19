require(data.table)
require(ggplot2)

args <- c("effstats.rds","fig_2.png")
args <- commandArgs(trailingOnly = TRUE)

stat.eff.dt <- readRDS(args[1])

vac.only <- stat.eff.dt[variable == "vac.eff", .(vac.eff=unique(med)), keyby=.(vac_mech, catchup, year)]
naive <- stat.eff.dt[variable == "ind.eff", .(assume.eff = med), keyby=.(vc_coverage, vac_mech, catchup, year)]

ribbon.dt <- naive[vac.only, on=.(vac_mech, catchup, year)]

combo.dt <- stat.eff.dt[variable == "combo.eff",.(eff = med), keyby=.(vc_coverage, vac_mech, catchup, year)]

other.ribbon <- combo.dt[naive, on=.(vc_coverage, vac_mech, catchup, year)][vac.only, on=.(vac_mech, catchup, year)]

other.ribbon[,
  lower := pmin(eff, assume.eff)
][,
  upper := pmax(eff, assume.eff)
][,
  outcome := ifelse(assume.eff < eff, "advantage", "disadvantage")
]

vac_mechs <- c("5th Sero", "Traditional")
vac_cols <- c("green", "blue")
names(vac_cols) <- vac_mechs
vacfact <- function(i) factor(vac_mechs[i+1], levels = rev(vac_mechs), ordered = T)

p <- ggplot(ribbon.dt) + aes(x=year, fill=factor(catchup)) +
  facet_grid(vacfact(vac_mech) ~ vc_coverage) +
  geom_ribbon(aes(ymin=vac.eff, ymax=assume.eff), alpha=0.4) +
  geom_line(data = combo.dt, mapping=aes(y=eff, color=factor(catchup), linetype="observed")) +
  geom_line(data = vac.only, mapping=aes(y=vac.eff, color=factor(catchup), linetype="vaccine only"), alpha=0.5) +
  geom_line(data = naive, mapping=aes(y=assume.eff, color=factor(catchup), linetype="naive est."), alpha=0.5) +
  theme_minimal() +
  scale_fill_manual(values=c(`0`="grey",`1`="black")) +
  scale_color_manual(values=c(`0`="darkgrey",`1`="black"))

# try drawing instead as area plot of obs around assumed
# when area above assumed, one color, another below
# crossing points need to be extracted

# generic problem: data series (x, yref, ycomp)
#  - need to insert crossing points
#  - crossing point defined as x1 -> x2 => (ycomp - yref) changes signs
#  - intersection occurs at some x, s.t.
#  (yref2-yref1)/(x2-x1)*x + bref = (ycomp2-ycomp1)/(x2-x1)*x + bcomp
#  bref = yref2 - (yref2-yref1)/(x2-x1)*x2, resp for bcomp
#  so x = (bcomp-bref)/(mref-mcomp)
#  - return list of geom_polygons

breakpoints <- other.ribbon[, {
  .(switches = cumsum(rle(outcome)$lengths))
}, by=.(vc_coverage, vac_mech, catchup)][
  switches != 40,
  .(leftyr=switches-1, rghtyr=switches),
  by=.(vc_coverage, vac_mech, catchup)
]

breakpoints[other.ribbon, 
  on=.(vc_coverage, vac_mech, catchup, leftyr==year)
]

other.ribbon[vc_coverage==25 & vac_mech == 0 & catchup == 1]

p2 <- ggplot(other.ribbon) +
  aes(x=year, fill=outcome, alpha = factor(vc_coverage)) +
  facet_grid(vacfact(vac_mech) ~ catchup) +
  geom_ribbon(aes(ymin=lower, ymax=upper)) +
  scale_fill_manual(values = c(advantage="green",disadvantage="red")) +
  scale_alpha_manual(values = c(`25`=.75,`50`=.50,`75`=.25)) +
  theme_minimal()

ggsave(
  tail(args,1), p, device = "png",
  width = 7.5, height = 5, dpi = "retina", units = "in"
)