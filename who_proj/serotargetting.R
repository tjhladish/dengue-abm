# seroprevalence targetting by beta
require(reshape2)
require(data.table)

readIn <- function(bmultipliers = seq(0.1, 2.0, len=20), src = "beta_dump-jan2016", refpop=15020) melt(
  setnames(fread(src), c("sample","beta_mul", paste0("year_",1:5))),
  id.vars = c("sample","beta_mul"),
  variable.name = "year", value.name = "seroprevalence"
)[,
  year := 94+as.integer(sub("year_","", year))
][,
  beta_times := bmultipliers[beta_mul+1]
][,
  seropositive := round(seroprevalence * refpop)
][,
  total := refpop
]

run_ref <- readIn()

plot(V1 ~ beta_times, data=run_ref[,list(mean(seropositive/total)), by=beta_times], ylim=c(0,1), xlim=c(0,3), type="b")
abline(h=c(.1,.3,.5,.7,.9), col="red")
abline(v=c(.48,.725,.975,1.375,2.8), col="blue")


# require(ggplot2)
# 
# ggplot(run_ref) + aes(group=sample, x=year, y=seroprevalence) + facet_grid(beta_mul ~ .) + geom_line(alpha=0.1)
# 
# ggplot(run_ref[, list(seroprevalence=mean(seroprevalence)), keyby=list(beta_mul)]) +
#   aes(x=bmultipliers[beta_mul+1], y=seroprevalence) + geom_point() +
#   stat_smooth(method = "glm", family=binomial(link="logit"))
# 
# glm.out <- glm(cbind(seropositive, total) ~ beta_times, data=run_ref[,list(seropositive=round(mean(seropositive))), keyby=list(beta_times,total)], family=binomial(logit))
# 
# seros <- c(.1, .3, .5, .7, .9)
# thing <- lm(sero~beta_mul, data=run_ref[, list(sero=mean(seroprevalence)), keyby=list(beta_mul)][,list(sero, beta_mul=beta_multiplier[beta_mul+1]^0.1)])
# newbeta <- ((seros - thing$coefficients[1])/thing$coefficients[2])^10
# 
# orig <- run_ref[, list(sero=mean(seroprevalence)), keyby=list(beta_mul)][, list(sero, beta=beta_multiplier[beta_mul+1], prov="original")]
# src <- rbind(orig, data.table(sero=seros, beta=newbeta, prov="fit"))
# ggplot(src) + aes(x=beta, y=sero, color=prov) + geom_point()
# 
# srcmod <- rbind(src, src[6:10][, list(sero, beta, prov="hand")])
# ggplot(srcmod) + aes(x=beta, y=sero, color=prov) + geom_point()
