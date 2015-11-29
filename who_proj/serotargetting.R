# seroprevalence targetting by beta
require(data.table)
require(reshape2)

beta_multiplier <- c(0.4, 0.8, 1, 1.2, 1.6)
run_ref <- melt(
  setnames(fread("beta_dump"), c("sample","beta_mul", paste0("year_",1:5))),
  id.vars = c("sample","beta_mul"),
  variable.name = "year", value.name = "seroprevalence"
)

run_ref[, year := 94+as.integer(sub("year_","", year))]
run_ref[, beta_times := beta_multiplier[beta_mul+1]]


require(ggplot2)

ggplot(run_ref) + aes(group=sample, x=year, y=seroprevalence) + facet_grid(beta_mul ~ .) + geom_line(alpha=0.1)

ggplot(run_ref[, mean(seroprevalence), keyby=list(beta_mul)]) + aes(x=beta_multiplier[beta_mul+1]^0.1, y=V1) + geom_point() + stat_smooth(method = "lm")

seros <- c(.1, .3, .5, .7, .9)
thing <- lm(sero~beta_mul, data=run_ref[, list(sero=mean(seroprevalence)), keyby=list(beta_mul)][,list(sero, beta_mul=beta_multiplier[beta_mul+1]^0.1)])
newbeta <- ((seros - thing$coefficients[1])/thing$coefficients[2])^10

orig <- run_ref[, list(sero=mean(seroprevalence)), keyby=list(beta_mul)][, list(sero, beta=beta_multiplier[beta_mul+1], prov="original")]
src <- rbind(orig, data.table(sero=seros, beta=newbeta, prov="fit"))
ggplot(src) + aes(x=beta, y=sero, color=prov) + geom_point()

srcmod <- rbind(src, src[6:10][, list(sero, beta, prov="hand")])
ggplot(srcmod) + aes(x=beta, y=sero, color=prov) + geom_point()
