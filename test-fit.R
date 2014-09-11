setwd("~/Documents/Research/age-specific-incidence/Spectrum/epp-spectrum-agefoi/")
load("R/south-africa-parameters.RData")
## source("R/spectrum.R")
source("C++/functions.R")
source("analysis-functions.R")
## source("R/incidence-model.R")
source("R/likelihood.R")

theta.start <- c(-27.2776051, 0.2286067, 10.5*c(1.5936094, 2.0196487 - 1.5936094, -0.4601538, -1.6417664, 1.4846658, 0.1725820, 0.3623600),
                 c(35, 20, 25, 15, 0.6))

param.start <- create.param(theta.start)
mod.start <- fnSpectrum(param.start)
round(prev(mod.start), 3)
cbind(round(prev(mod.start, a15to49.idx, m.idx), 3),
      round(prev(mod.start, a15to49.idx, f.idx), 3))

ll(theta.start)
anc.dat$logit.var <- sa.anc$logit.var * 100
ll(theta.start)

test.fit.infl <- optim(theta.start, ll, DEBUG=TRUE, method="Nelder-Mead", control=list(fnscale=-1, maxit=10000, trace=4, REPORT=1))
test.fit.infl.b <- optim(test.fit.infl$par, ll, DEBUG=TRUE, method="Nelder-Mead", control=list(fnscale=-1, maxit=10000, trace=4, REPORT=1))

anc.dat$logit.var <- sa.anc$logit.var
ll(test.fit.infl.b$par)
test.fit <- optim(test.fit.infl.b$par, ll, DEBUG=TRUE, method="Nelder-Mead", control=list(fnscale=-1, maxit=10000, trace=4, REPORT=1))
test.fit.b <- optim(test.fit$par, ll, DEBUG=TRUE, method="Nelder-Mead", control=list(fnscale=-1, maxit=10000, trace=4, REPORT=1))
test.fit.c <- optim(test.fit.b$par, ll, DEBUG=TRUE, method="Nelder-Mead", control=list(fnscale=-1, maxit=10000, trace=4, REPORT=1))

save(test.fit.infl, test.fit.infl.b, test.fit, test.fit.b, test.fit.c, file="test-fit-workspace.RData")

########################
####  Plot results  ####
########################

incid <- function(mod, param, age.idx=a15to49.idx, sex.idx=c(m.idx, f.idx)){
  rval <- param$rVec[match(1:43+0.5, proj.steps)]
  return(sapply(1:43, function(i){ age.inc <- fnAgeInc(mod[i,,,,], rval[i], 0, param)
                                   return(sum(mod[i,sex.idx,age.idx,1,1] * age.inc[sex.idx,age.idx]) / sum(mod[i,sex.idx,age.idx,1,1])) }))
}

param.fit <- create.param(test.fit.c$par)
mod.fit <- fnSpectrum(param.fit)

sa.ageprev$years <- c(2002, 2005, 2008, 2012)

quartz(w=5, h=2, pointsize=9)

par(mfrow=c(1, 2), mar=c(2, 3.5, 2, 3.5), cex=0.8, mgp=c(2, 0.5, 0), tcl=-0.25)
plot(1970:2012, prev(mod.fit), type="l", ylim=c(0, 0.25), xlim=c(1990, 2012),
     ylab="Prevalence, 15-49 y")
points(sa.ageprev$years, sa.15to49$prevalence)
segments(sa.ageprev$years, y0=sa.15to49$lowerCI, y1=sa.15to49$upperCI)
par(new=TRUE)
plot(1970:2012, incid(mod.fit, param.fit), type="l", ylim=c(0, 0.035), xlim=c(1990, 2012), axes=FALSE, lty=2, ylab="", xlab="", bty="n")
axis(4)
mtext("Incidence, 15-49 y", 4, 2, cex=0.8)
####
matplot(1970:2012, cbind(prev(mod.fit, sex=m.idx), prev(mod.fit, sex=f.idx)), type="l", lty=1, ylim=c(0, 0.25), xlim=c(1990, 2012), col=3:4,
        ylab="Prevalence, 15-49 y")
matpoints(sa.ageprev$years, cbind(sa.male15to49$prevalence, sa.female15to49$prevalence), pch=1, col=3:4)
segments(sa.ageprev$years, y0=sa.male15to49$lowerCI, y1=sa.male15to49$upperCI, col=3)
segments(sa.ageprev$years, y0=sa.female15to49$lowerCI, y1=sa.female15to49$upperCI, col=4)
par(new=TRUE)
matplot(1970:2012, cbind(incid(mod.fit, param.fit, sex=m.idx), incid(mod.fit, param.fit, sex=f.idx)), type="l", lty=2, ylim=c(0, 0.04), xlim=c(1990, 2012), col=3:4,
        ylab="", axes=FALSE, bty="n")
axis(4)
mtext("Incidence, 15-49 y", 4, 2, cex=0.8)
        

####  Age-specific prevalence

ageprev.fit <- ageprev(mod.fit[hhsurv.ageprev$idx,,,,], age15plus.idx, c(1:9, rep(10, 5)))

quartz(w=5, h=7, pointsize=9)

par(mfrow=c(4, 2), mar=c(2, 3.5, 1, 1), mgp=c(2, 0.5, 0), tcl=-0.25, cex=0.8)

for(surv in 1:4)
  for(sex in 1:2){
    plot(3:12*5, ageprev.fit[surv,sex,], type="l", ylim=c(0, 0.4), lty=1, lwd=2, xlab="", ylab="Prevalence")
    points(3:12*5, t(sa.ageprev$prev[surv,sex,]), col=1, pch=1, type="b", lty=2)
    segments(3:12*5, y0=sa.ageprev$lowerCI[surv,sex,], y1=sa.ageprev$upperCI[surv,sex,], col=1, lty=2)
    mtext(paste(c("Men", "Women")[sex], "-", c("2002", "2005", "2008", "2012")[surv]), 3, -1.5, adj=0.95, font=2, cex=0.8)
}
