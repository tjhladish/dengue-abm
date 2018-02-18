rm(list = ls())
require(deSolve)
require(data.table)
require(ggplot2)
require(gridExtra)

# define three SIR models, ideally in terms of each other

SIR <- function(t, SIRt, params) with(params, {
  S <- SIRt[1]; I <- SIRt[2]; R <- SIRt[3]
  infection <- b*S*I
  recovery <- lambda*I
  dS <- -infection
  dI <- infection - recovery
  dR <- recovery
  return(list(c(dS, dI, dR)))
})

waningSIR <- function(t, SIRt, params) {
  dSIR <- SIR(t, SIRt, params)
  R <- SIRt[3]
  waning <- params$tau*R
  dSIR[[1]][1] <- dSIR[[1]][1] + waning
  dSIR[[1]][3] <- dSIR[[1]][3] - waning
  return(dSIR)
}

seasSIR <- function(t, SIRt, params) {
  params$b <- params$b*params$sf(t)
  return(waningSIR(t, SIRt, params))
}

## given functions:
## - run 4 sims:
##  * 20 for baseline
##  * 15 starting from baseline 5 for 3x levels

runbaseline <- function(model, pars, nm="baseline", SIR0=c(S=0.99, I=0.01, R=0)) {
  res <- data.table(ode(SIR0, seq(0, 20, by=0.001), model, pars))
  res[, name := nm ]
}

runintervention <- function(baseline, model, pars, nm, start = 3) {
  SIR0 <- as.numeric(baseline[time == start,c("S","I","R"),with=F])
  names(SIR0) <- c("S","I","R")
  res <- runbaseline(model, pars, nm, SIR0)
  res[, time := time + start ]
  res
}

runprog <- function(model, pars, sq=((1:3)/4)*100, start = 3) {
  basebeta <- pars$b
  bs <- runbaseline(model, pars)
  res <- bs
  for (mod in sq) {
    pars$b <- basebeta * (1-(mod/100))
    add <- runintervention(bs, model, pars, mod, start)
    # browser()
    res <- rbind(res, add)
  }
  res <- melt(res, id.vars = c("name","time"))
  res[, name := factor(name, levels = c("baseline",sq), ordered = T)]
}

lamb <- 52 # 1 per week == 52 per year
Beta <- 10/1*52 # 10 new cases per infectious person over one week = 10 per (2/52s year)

SIRpars <- list(b=Beta, lambda=lamb)

SIR.dt <- runprog(model=SIR, pars=SIRpars, start=0)

ggplot(SIR.dt) + aes(x=time, linetype=name, color=variable, y=value) + geom_line() +
  coord_cartesian(xlim=c(0,0.25))

tau <- 1/2

waningSIRpars <- list(b=Beta, lambda=lamb, tau=tau)
waningSIR.dt <- runprog(model=waningSIR, pars=waningSIRpars)

ggplot(waningSIR.dt) + aes(x=time, linetype=name, color=variable, y=value) + geom_line() +
  coord_cartesian(xlim=c(0,10))

sigmoid <- function(x,L=1) L/(1+exp(-10*(sin(2*pi*x))))

seasSIRpars <- list(b=Beta, lambda=lamb, tau=tau, sf=sigmoid)
seasSIR.dt <- runprog(model=seasSIR, pars=seasSIRpars)

ggplot(seasSIR.dt) + aes(x=time, linetype=name, color=variable, y=value) + geom_line() +
  coord_cartesian(xlim=c(0,10))

ggplot_res <- function(dt, tlim) {
  slice <- dt[time <= tlim]
  lvls <- slice[, unique(name)]
  p1 <- ggplot(slice) + aes(x=time, linetype=name, color=variable, y=value) + geom_line()
  p2 <- ggplot(wide) + aes(x=I, color=name, y=R) + geom_path() + scale_color_discrete(breaks=lvls, drop=F)
  wide <- dcast(slice, time + name ~ variable)
  base <- wide[name == "baseline"]
  iv <- wide[name != "baseline"]
  
  eff.dt <- iv[base, on="time"][,.(
    base.I = i.I, I,
    cum.base=cumsum(i.I), cum.I=cumsum(I),
    time
    ), by="name"
  ][, 
    .(eff=(base.I-I)/base.I, cum.eff=(cum.base-cum.I)/cum.base),
    keyby=.(name,time)
  ]
  p3 <- ggplot(eff.dt) + aes(x=time, color=name, y=eff) + geom_line() + scale_color_discrete(breaks=lvls[-1], drop=F)
  p4 <- ggplot(eff.dt) + aes(x=time, color=name, y=cum.eff) + geom_line() + scale_color_discrete(breaks=lvls[-1], drop=F)
  grid.arrange(p1, p2, p3, p4, nrow = 2)
}

ggplot_res(SIR.dt, 0.25)
ggplot_res(waningSIR.dt, 10)
ggplot_res(seasSIR.dt, 20)

# introduction - once per month, 1% of pop successful exposed to dengue?
# so dSintro <- -sigma*S, sigma = 0.01*12

SIR <- function(t, SIR, params) with(params, {
  S <- SIR[1]; I <- SIR[2]; R <- SIR[3]
  infection <- b*S*I
  recovery <- lambda*I
  dS <- -infection
  dI <- infection - recovery
  dR <- recovery
  return(list(c(dS, dI, dR)))
})

tauSIR <- function(t, SIR, params) with(params, {
  S <- SIR[1]; I <- SIR[2]; R <- SIR[3]
  infection <- b*S*I
  recovery <- lambda*I
  waning <- tau*R
  dS <- -infection + waning
  dI <- infection - recovery
  dR <- recovery - waning
  return(list(c(dS, dI, dR)))
})

seasTauSIR <- function(t, SIR, params) with(params, {
  S <- SIR[1]; I <- SIR[2]; R <- SIR[3];
  infection <- b*sf(t)*S*I
  recovery <- lambda*I
  waning <- tau*R
  dS <- -infection + waning
  dI <- infection - recovery
  dR <- recovery - waning
  return(list(c(dS, dI, dR)))
})

sigmoid <- function(x,L=1) L/(1+exp(-10*(sin(2*pi*x))))
plot(seq(-2,2,by=0.01), sigmoid(seq(-2,2,by=0.01)))

lamb <- 52 # 1 per week == 52 per year
Beta <- 10/1*52 # 10 new cases per infectious person over one week = 10 per (2/52s year)

genfunc <- function(f,p,end=5,SIR0=c(S=0.99,I=0.01,R=0)) melt(
  data.table(ode(SIR0, seq(0,end,by=0.001), f, p)),
  id.vars = "time"
)

ddegenfunc <- function(f,p,switchF,mapF,end=5) melt(
  data.table(dde(
    c(S=0.99,I=0.01,R=0), seq(0,end,by=0.001), f, p
  )), id.vars = "time"
)

SIR.dt <- genfunc(SIR, list(b=Beta,lambda=lamb))

ggplot(SIR.dt) + aes(x=time, y=value, color=variable) + geom_line()

tau <- 1/2

tauSIR.dt <- genfunc(tauSIR, list(b=Beta,lambda=lamb,tau=tau),10)

ggplot(tauSIR.dt) + aes(x=time, y=value, color=variable) + geom_line()

seasTauSIR.dt <- genfunc(seasTauSIR, list(b=Beta,lambda=lamb,tau=tau,sf=sigmoid),20)

ggplot(seasTauSIR.dt) + aes(x=time, y=value, color=variable) + geom_line() +
  scale_x_continuous(breaks = 0:20)

wide <- dcast(seasTauSIR.dt, time ~ variable)

ggplot(wide) + aes(x=I, y=R) + geom_path(alpha=0.5)

seasTauSIR2.dt <- genfunc(seasTauSIR, list(b=Beta,lambda=lamb,tau=tau,sf=sigmoid),5)
midSIR0 <- seasTauSIR2.dt[time == 5,{res<-as.vector(value); names(res) <- as.vector(variable); res }]

seasTauSIR0.dt <- genfunc(seasTauSIR, list(b=Beta*0.2,lambda=lamb,tau=tau,sf=sigmoid),15,midSIR0)
seasTauSIR3.dt <- genfunc(seasTauSIR, list(b=Beta*0.3,lambda=lamb,tau=tau,sf=sigmoid),15,midSIR0)
seasTauSIR4.dt <- genfunc(seasTauSIR, list(b=Beta*0.4,lambda=lamb,tau=tau,sf=sigmoid),15,midSIR0)
seasTauSIR5.dt <- genfunc(seasTauSIR, list(b=Beta*0.5,lambda=lamb,tau=tau,sf=sigmoid),15,midSIR0)
seasTauSIR6.dt <- genfunc(seasTauSIR, list(b=Beta*0.6,lambda=lamb,tau=tau,sf=sigmoid),15,midSIR0)

stitch <- rbind(seasTauSIR2.dt[, level := "pre"],
  seasTauSIR3.dt[time!=0,list(time=time+5, variable, value, level="post-70")],
  seasTauSIR4.dt[time!=0,list(time=time+5, variable, value, level="post-60")],
  seasTauSIR5.dt[time!=0,list(time=time+5, variable, value, level="post-50")]
)

wide <- dcast(stitch, level + time ~ variable)

ggplot(wide) + aes(x=I, y=R, color=level) + geom_path(alpha=0.5)

ggplot(stitch) + aes(x=time, y=value, color=variable) + geom_line() +
  scale_x_continuous(breaks = 0:20)
