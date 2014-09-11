setwd("~/Documents/Research/age-specific-incidence/Spectrum/epp-spectrum-agefoi/")
load("R/south-africa-parameters.RData")
source("C++/functions.R")
source("analysis-functions.R")
source("R/likelihood.R")
source("R/incidence-model.R")

source("R/IMIS.R")

anc.dat$logit.var <- sa.anc$logit.var * 100

set.seed(464417)
test.imis1 <- IMIS(2e7, 10000, 3000, 500, 0)

set.seed(5673167)
test.imis2 <- IMIS(2e7, 10000, 3000, 500)

set.seed(37831476)
test.imis3 <- IMIS(2e7, 10000, 3000, 500)

set.seed(55556272)
test.imis4 <- IMIS(2e7, 10000, 3000, 500)

save(test.imis1, test.imis2, test.imis3, test.imis4,
     file="test-imis.RData")

MODEL <- "rmat"
set.seed(58559308); test.rmat1 <- IMIS(2e6, 10000, 3000, 500)
set.seed(10389276); test.rmat2 <- IMIS(2e6, 10000, 3000, 500)
set.seed(19403255); test.rmat3 <- IMIS(2e6, 10000, 3000, 500)
set.seed(89656651); test.rmat4 <- IMIS(2e6, 10000, 3000, 500)


MODEL <- "incrr"
set.seed(8993068); test.incrr1 <- IMIS(2e6, 10000, 3000, 500)
set.seed(29733230); test.incrr2 <- IMIS(2e6, 10000, 3000, 500)
set.seed(75238206); test.incrr3 <- IMIS(2e6, 10000, 3000, 500)
set.seed(2588340); test.incrr4 <- IMIS(2e6, 10000, 3000, 500)

save(test.rmat1, test.rmat2, test.rmat3, test.rmat4, 
     test.incrr1, test.incrr2, test.incrr3, test.incrr4,
     file="test-imis-rmat-incrr.RData")


############################################################
############################################################

unaids2013.sa.prev <- data.frame(mean=c(0.3, 0.6, 1.1, 2.0, 3.3, 5.0, 7.0, 9.2, 11.3, 13.2, 14.9, 16.2, 17.1, 17.8, 18.3, 18.5, 18.7, 18.8, 18.9, 18.9, 18.9, 18.8, 18.9, 19.1),
                                 lowerCI=c(0.3, 0.5, 1.0, 1.9, 3.1, 4.7, 6.6, 8.7, 10.8, 12.6, 14.2, 15.5, 16.4, 17.1, 17.5, 17.8, 17.9, 18.0, 18.1, 18.1, 18.0, 18.0, 18.0, 18.1),
                                 upperCI=c(0.3, 0.6, 1.2, 2.2, 3.5, 5.3, 7.4, 9.7, 11.9, 13.8, 15.5, 16.8, 17.9, 18.6, 19.0, 19.3, 19.5, 19.6, 19.7, 19.8, 19.7, 19.7, 19.8, 19.9),
                                 row.names=1990:2013)/100

unaids2013.sa.incid <- data.frame(mean=c(0.17, 0.33, 0.59, 0.98, 1.46, 1.99, 2.48, 2.83, 3.01, 3.04, 2.95, 2.78, 2.59, 2.40, 2.24, 2.11, 2.01, 1.93, 1.88, 1.84, 1.79, 1.69, 1.52, 1.36),
                                lowerCI=c(0.15, 0.30, 0.55, 0.92, 1.38, 1.88, 2.36, 2.70, 2.89, 2.92, 2.83, 2.68, 2.49, 2.31, 2.15, 2.02, 1.92, 1.85, 1.79, 1.76, 1.70, 1.60, 1.43, 1.26),
                                upperCI=c(0.18, 0.36, 0.64, 1.05, 1.56, 2.11, 2.61, 2.96, 3.13, 3.15, 3.05, 2.88, 2.68, 2.50, 2.33, 2.20, 2.10, 2.02, 1.96, 1.93, 1.88, 1.79, 1.62, 1.45),
                                row.names=1990:2013)/100

unaids2013.sa.f15to24prev <- data.frame(mean=c(0.4, 0.8, 1.6, 2.8, 4.7, 7.0, 9.7, 12.5, 14.9, 16.9, 18.3, 19.1, 19.4, 19.2, 18.7, 18.1, 17.3, 16.6, 15.9, 15.3, 14.8, 14.3, 13.8, 13.1),
                                        lowerCI=c(0.3, 0.7, 1.4, 2.5, 4.1, 6.2, 8.7, 11.3, 13.6, 15.4, 16.7, 17.4, 17.6, 17.5, 17.0, 16.4, 15.7, 15.0, 14.4, 13.9, 13.4, 13.0, 12.5, 11.9),
                                        upperCI=c(0.5, 1.0, 1.8, 3.3, 5.4, 8.1, 11.2, 14.5, 17.4, 19.9, 21.7, 22.8, 23.2, 23.2, 22.8, 22.0, 21.1, 20.3, 19.4, 18.7, 18.0, 17.4, 16.8, 16.1)
                                        row.names=1990:2013)

unaids2013.sa.m15to24prev <- data.frame(mean=c(0.2, 0.3, 0.6, 1.1, 1.7, 2.5, 3.4, 4.2, 4.9, 5.4, 5.8, 5.9, 5.9, 5.7, 5.5, 5.2, 5.0, 4.8, 4.6, 4.5, 4.4, 4.2, 4.1, 4.0),
                                        lowerCI=c(0.1, 0.2, 0.4, 0.7, 1.2, 1.7, 2.3, 2.8, 3.2, 3.6, 3.7, 3.7, 3.7, 3.6, 3.4, 3.2, 3.1, 3.0, 2.9, 2.8, 2.7, 2.7, 2.6, 2.5),
                                        upperCI=c(0.2, 0.4, 0.8, 1.3, 2.1, 3.1, 4.2, 5.3, 6.3, 7.1, 7.8, 8.2, 8.5, 8.5, 8.3, 8.0, 7.7, 7.4, 7.1, 6.8, 6.6, 6.4, 6.1, 5.9),
                                        row.names=1990:2013)

############################################################
############################################################

incid <- function(mod, param, age.idx=a15to49.idx, sex.idx=c(m.idx, f.idx)){
  rval <- param$rVec[match(1:43+0.5, proj.steps)]
  return(sapply(1:43, function(i){ age.inc <- fnAgeInc(mod[i,,,,], rval[i], 0, param)
                                   return(sum(mod[i,sex.idx,age.idx,1,1] * age.inc[sex.idx,age.idx]) / sum(mod[i,sex.idx,age.idx,1,1])) }))
}

############################################################
############################################################

load("test-imis-rmat-incrr.RData")

MODEL <- "rmat"
round(prev(fnSpectrum(create.param(test.rmat1$resample[8,]))), 3)
round(prev(fnSpectrum(create.param(test.rmat2$resample[8,]))), 3)
round(prev(fnSpectrum(create.param(test.rmat3$resample[8,]))), 3)
round(prev(fnSpectrum(create.param(test.rmat4$resample[1,]))), 3)

param.rmat1 <- apply(test.rmat1$resample, 1, create.param)
param.rmat2 <- apply(test.rmat2$resample, 1, create.param)
param.rmat3 <- apply(test.rmat3$resample, 1, create.param)
param.rmat4 <- apply(test.rmat4$resample, 1, create.param)

mod.rmat1 <- lapply(param.rmat1, fnSpectrum)
mod.rmat2 <- lapply(param.rmat2, fnSpectrum)
mod.rmat3 <- lapply(param.rmat3, fnSpectrum)
mod.rmat4 <- lapply(param.rmat4, fnSpectrum)

ll.rmat1 <- apply(test.rmat1$resample, 1, ll)
ll.rmat2 <- apply(test.rmat2$resample, 1, ll)
ll.rmat3 <- apply(test.rmat3$resample, 1, ll)
ll.rmat4 <- apply(test.rmat4$resample, 1, ll)

lprior.rmat1 <- apply(test.rmat1$resample, 1, lprior)
lprior.rmat2 <- apply(test.rmat2$resample, 1, lprior)
lprior.rmat3 <- apply(test.rmat3$resample, 1, lprior)
lprior.rmat4 <- apply(test.rmat4$resample, 1, lprior)



MODEL <- "incrr"
param.incrr1 <- apply(test.incrr1$resample, 1, create.param)
param.incrr2 <- apply(test.incrr2$resample, 1, create.param)
param.incrr3 <- apply(test.incrr3$resample, 1, create.param)
param.incrr4 <- apply(test.incrr4$resample, 1, create.param)

mod.incrr1 <- lapply(param.incrr1, fnSpectrum)
mod.incrr2 <- lapply(param.incrr2, fnSpectrum)
mod.incrr3 <- lapply(param.incrr3, fnSpectrum)
mod.incrr4 <- lapply(param.incrr4, fnSpectrum)

ll.incrr1 <- apply(test.incrr1$resample, 1, ll)
ll.incrr2 <- apply(test.incrr2$resample, 1, ll)
ll.incrr3 <- apply(test.incrr3$resample, 1, ll)
ll.incrr4 <- apply(test.incrr4$resample, 1, ll)

lprior.incrr1 <- apply(test.incrr1$resample, 1, lprior)
lprior.incrr2 <- apply(test.incrr2$resample, 1, lprior)
lprior.incrr3 <- apply(test.incrr3$resample, 1, lprior)
lprior.incrr4 <- apply(test.incrr4$resample, 1, lprior)


########################
####  Plot results  ####
########################

incid <- function(mod, param, age.idx=a15to49.idx, sex.idx=c(m.idx, f.idx)){
  rval <- param$rVec[match(1:43+0.5, proj.steps)]
  return(sapply(1:43, function(i){ age.inc <- fnAgeInc(mod[i,,,,], rval[i], 0, param)
                                   return(sum(mod[i,sex.idx,age.idx,1,1] * age.inc[sex.idx,age.idx]) / sum(mod[i,sex.idx,age.idx,1,1])) }))
}


prev.rmat1 <- sapply(mod.rmat1, prev)
prev.rmat2 <- sapply(mod.rmat2, prev)
prev.rmat3 <- sapply(mod.rmat3, prev)
prev.rmat4 <- sapply(mod.rmat4, prev)

prev.incrr1 <- sapply(mod.incrr1, prev)
prev.incrr2 <- sapply(mod.incrr2, prev)
prev.incrr3 <- sapply(mod.incrr3, prev)
prev.incrr4 <- sapply(mod.incrr4, prev)

matplot(1970:2012, sapply(mod.rmat1, prev))

incid.rmat1 <- mapply(incid, mod.rmat1, param.rmat1)
incid.rmat2 <- mapply(incid, mod.rmat2, param.rmat2)
incid.rmat3 <- mapply(incid, mod.rmat3, param.rmat3)
incid.rmat4 <- mapply(incid, mod.rmat4, param.rmat4)

incid.incrr1 <- mapply(incid, mod.incrr1, param.incrr1)
incid.incrr2 <- mapply(incid, mod.incrr2, param.incrr2)
incid.incrr3 <- mapply(incid, mod.incrr3, param.incrr3)
incid.incrr4 <- mapply(incid, mod.incrr4, param.incrr4)


save("incid.rmat1", "incid.rmat2", "incid.rmat3", "incid.rmat4",
     "ll.rmat1",  "ll.rmat2",  "ll.rmat3",  "ll.rmat4",
     "lprior.rmat1", "lprior.rmat2", "lprior.rmat3", "lprior.rmat4",
     ## "mod.rmat1", "mod.rmat2", "mod.rmat3", "mod.rmat4",
     ## "param.rmat1", "param.rmat2", "param.rmat3", "param.rmat4",
     "prev.rmat1", "prev.rmat2", "prev.rmat3", "prev.rmat4",
     "test.rmat1", "test.rmat2", "test.rmat3", "test.rmat4",
     "incid.incrr1", "incid.incrr2", "incid.incrr3", "incid.incrr4",
     "ll.incrr1",  "ll.incrr2",  "ll.incrr3",  "ll.incrr4",
     "lprior.incrr1", "lprior.incrr2", "lprior.incrr3", "lprior.incrr4",
     ## "mod.incrr1", "mod.incrr2", "mod.incrr3", "mod.incrr4",
     ## "param.incrr1", "param.incrr2", "param.incrr3", "param.incrr4",
     "prev.incrr1", "prev.incrr2", "prev.incrr3", "prev.incrr4",
     "test.incrr1", "test.incrr2", "test.incrr3", "test.incrr4",
     file="test-imis-analysis.RData")


############################################################
############################################################

load("test-imis-analysis.RData")

sa.ageprev$years <- c(2002, 2005, 2008, 2012)


quartz(w=8, h=4, pointsize=10)
par(mfrow=c(1, 2), mar=c(2, 2.5, 2.5, 0.5), mgp=c(2, 0.5, 0), tcl=-0.25, las=1)

plot(1970:2012, rowMeans(prev.rmat1),
        lty=1, lwd=1.5, col=4, type="l", xlim=c(1990, 2012), ylim=c(0, 0.2),
     main="FOI matrix", ylab="", xlab="")
matlines(1970:2012, t(apply(prev.rmat1, 1, quantile, c(0.025, 0.975))), 
     lty=2, lwd=1.5, col=4)
points(sa.ageprev$years, sa.15to49$prevalence, lwd=1.5)
segments(sa.ageprev$years, y0=sa.15to49$lowerCI, y1=sa.15to49$upperCI, lwd=1.5)
##
plot(1970:2012, rowMeans(prev.incrr1),
        lty=1, lwd=1.5, col=2, type="l", xlim=c(1990, 2012), ylim=c(0, 0.2),
     main="Incidence RR", ylab="", xlab="")
matlines(1970:2012, t(apply(prev.incrr1, 1, quantile, c(0.025, 0.975))), 
     lty=2, lwd=1.5, col=2)
points(sa.ageprev$years, sa.15to49$prevalence, lwd=1.5)
segments(sa.ageprev$years, y0=sa.15to49$lowerCI, y1=sa.15to49$upperCI, lwd=1.5)


###################################
####  HIV incidence over time  ####
###################################

quartz(w=8, h=4, pointsize=10)
par(mfrow=c(1, 2), mar=c(2, 2.5, 2.5, 0.5), mgp=c(2, 0.5, 0), tcl=-0.25, las=1)

plot(1970:2012, rowMeans(incid.rmat1),
     lty=1, lwd=1.5, col=4, type="l", xlim=c(1990, 2012), ylim=c(0, 0.035),
     main="FOI matrix", ylab="", xlab="")
matlines(1970:2012, t(apply(incid.rmat1, 1, quantile, c(0.025, 0.975))), 
     lty=2, lwd=1.5, col=4)
##
plot(1970:2012, rowMeans(incid.incrr1),
        lty=1, lwd=1.5, col=2, type="l", xlim=c(1990, 2012), ylim=c(0, 0.035),
     main="Incidence RR", ylab="", xlab="")
matlines(1970:2012, t(apply(incid.incrr1, 1, quantile, c(0.025, 0.975))), 
     lty=2, lwd=1.5, col=2)


#############################
####  Prevalence by sex  ####
#############################

maleprev.rmat1 <- sapply(mod.rmat1, prev, sex=m.idx)
femaleprev.rmat1 <- sapply(mod.rmat1, prev, sex=f.idx)
maleprev.incrr1 <- sapply(mod.incrr1, prev, sex=m.idx)
femaleprev.incrr1 <- sapply(mod.incrr1, prev, sex=f.idx)

quartz(w=8, h=4, pointsize=10)

par(mfrow=c(1, 2), mar=c(2, 2.5, 2.5, 0.5), mgp=c(2, 0.5, 0), tcl=-0.25, las=1)
matplot(1970:2012, cbind(rowMeans(maleprev.rmat1), rowMeans(femaleprev.rmat1)), type="l", lty=1, ylim=c(0, 0.25), xlim=c(1990, 2012), col=3:4, lwd=1.5,
        main="FOI matrix", ylab="", xlab="")
matlines(1970:2012, t(rbind(apply(maleprev.rmat1, 1, quantile, c(0.025, 0.975)), apply(femaleprev.rmat1, 1, quantile, c(0.025, 0.975)))), lty=2, lwd=1.5, col=c(3,3,4,4))
matpoints(sa.ageprev$years, cbind(sa.male15to49$prevalence, sa.female15to49$prevalence), pch=1, col=3:4, lwd=1.5)
segments(sa.ageprev$years, y0=sa.male15to49$lowerCI, y1=sa.male15to49$upperCI, col=3, lwd=1.5)
segments(sa.ageprev$years, y0=sa.female15to49$lowerCI, y1=sa.female15to49$upperCI, col=4, lwd=1.5)
####
matplot(1970:2012, cbind(rowMeans(maleprev.incrr1), rowMeans(femaleprev.incrr1)), type="l", lty=1, ylim=c(0, 0.25), xlim=c(1990, 2012), col=3:4, lwd=1.5,
        main="Incidence RR", ylab="", xlab="")
matlines(1970:2012, t(rbind(apply(maleprev.incrr1, 1, quantile, c(0.025, 0.975)), apply(femaleprev.incrr1, 1, quantile, c(0.025, 0.975)))), lty=2, lwd=1.5, col=c(3,3,4,4))
matpoints(sa.ageprev$years, cbind(sa.male15to49$prevalence, sa.female15to49$prevalence), pch=1, col=3:4, lwd=1.5)
segments(sa.ageprev$years, y0=sa.male15to49$lowerCI, y1=sa.male15to49$upperCI, col=3, lwd=1.5)
segments(sa.ageprev$years, y0=sa.female15to49$lowerCI, y1=sa.female15to49$upperCI, col=4, lwd=1.5)

###################################
####  Age specific prevalence  ####
###################################


sim <- mod.incrr1[[1]]


ageprev.incrr1 <- mclapply(mod.incrr1, function(sim) ageprev(sim[hhsurv.ageprev$idx,,,,], age15plus.idx, c(1:9, rep(10, 5))))
ageprev.rmat1 <- mclapply(mod.rmat1, function(sim) ageprev(sim[hhsurv.ageprev$idx,,,,], age15plus.idx, c(1:9, rep(10, 5))))

  


quartz(w=8, h=6.5, pointsize=9)

par(mfrow=c(4, 4), mar=c(2, 2, 0.5, 0.5), mgp=c(2, 0.5, 0), tcl=-0.25, cex=0.8, las=1)
##
for(sex in 1:2)
for(surv in 1:4){
    plot(3:12*5, rowMeans(sapply(ageprev.rmat1, "[", i=surv, j=sex, TRUE)), type="l", ylim=c(0, 0.4), lty=1, lwd=1.5, xlab="", ylab="", col=4)
    matlines(3:12*5, t(apply(sapply(ageprev.rmat1, "[", i=surv, j=sex, TRUE), 1, quantile, c(0.025, 0.975))), lty=2, lwd=1.5, col=4)
    points(3:12*5, t(sa.ageprev$prev[surv,sex,]), col=1, pch=1, lwd=1.5)
    segments(3:12*5, y0=sa.ageprev$lowerCI[surv,sex,], y1=sa.ageprev$upperCI[surv,sex,], col=1, lty=1, lwd=1.5)
    mtext(paste(c("Men", "Women")[sex], "-", c("2002", "2005", "2008", "2012")[surv]), 3, -1.5, adj=0.95, font=2, cex=0.8)
}
##
for(sex in 1:2)
for(surv in 1:4){
    plot(3:12*5, rowMeans(sapply(ageprev.incrr1, "[", i=surv, j=sex, TRUE)), type="l", ylim=c(0, 0.4), lty=1, lwd=1.5, xlab="", ylab="", col=2)
    matlines(3:12*5, t(apply(sapply(ageprev.incrr1, "[", i=surv, j=sex, TRUE), 1, quantile, c(0.025, 0.975))), lty=2, lwd=1.5, col=2)
    points(3:12*5, t(sa.ageprev$prev[surv,sex,]), col=1, pch=1, lwd=1.5)
    segments(3:12*5, y0=sa.ageprev$lowerCI[surv,sex,], y1=sa.ageprev$upperCI[surv,sex,], col=1, lty=1, lwd=1.5)
    mtext(paste(c("Men", "Women")[sex], "-", c("2002", "2005", "2008", "2012")[surv]), 3, -1.5, adj=0.95, font=2, cex=0.8)
}


##########################################
####  Age specific prevalence trends  ####
##########################################

load("~/Google Drive/Modelling Consortium/Work Package 12 (Model validation)/HIVMC12.2 sa-hsrc-comparison/hsrc data/hsrc-estimates-raw-data.RData")


a25to49.idx <- 6:10
a50plus.idx <- 11:17

male15to24prev.rmat1 <- sapply(mod.rmat1, prev, age=a15to24.idx, sex=m.idx)
male25to49prev.rmat1 <- sapply(mod.rmat1, prev, age=a25to49.idx, sex=m.idx)
male50plusprev.rmat1 <- sapply(mod.rmat1, prev, age=a50plus.idx, sex=m.idx)

female15to24prev.rmat1 <- sapply(mod.rmat1, prev, age=a15to24.idx, sex=f.idx)
female25to49prev.rmat1 <- sapply(mod.rmat1, prev, age=a25to49.idx, sex=f.idx)
female50plusprev.rmat1 <- sapply(mod.rmat1, prev, age=a50plus.idx, sex=f.idx)

male15to24prev.incrr1 <- sapply(mod.incrr1, prev, age=a15to24.idx, sex=m.idx)
male25to49prev.incrr1 <- sapply(mod.incrr1, prev, age=a25to49.idx, sex=m.idx)
male50plusprev.incrr1 <- sapply(mod.incrr1, prev, age=a50plus.idx, sex=m.idx)

female15to24prev.incrr1 <- sapply(mod.incrr1, prev, age=a15to24.idx, sex=f.idx)
female25to49prev.incrr1 <- sapply(mod.incrr1, prev, age=a25to49.idx, sex=f.idx)
female50plusprev.incrr1 <- sapply(mod.incrr1, prev, age=a50plus.idx, sex=f.idx)

save("male15to24prev.rmat1", "male25to49prev.rmat1", "male50plusprev.rmat1",
     "female15to24prev.rmat1", "female25to49prev.rmat1", "female50plusprev.rmat1",
     "male15to24prev.incrr1", "male25to49prev.incrr1", "male50plusprev.incrr1",
     "female15to24prev.incrr1", "female25to49prev.incrr1", "female50plusprev.incrr1",
     file="test-imis-age-specific-trends.RData")


quartz(h=4.5, w=8, pointsize=10)

par(mfrow=c(2, 3), mar=c(2, 2.5, 0.5, 0.5), mgp=c(2, 0.5, 0), tcl=-0.25, cex=1, las=1)
##
plot(1970:2012, rowMeans(male15to24prev.rmat1), type="l", lty=1, lwd=1.5, col=4, xlim=c(1995, 2012), ylim=c(0, 0.15), ylab="", xlab="")
matlines(1970:2012, t(apply(male15to24prev.rmat1, 1, quantile, c(0.025, 0.975))), lty=2, lwd=1.5, col=4)
points(sa.ageprev$years, hsrc.male15to24.prev$mean, lwd=1.5, cex=0.8)
segments(sa.ageprev$years, y0=hsrc.male15to24.prev$lowerCI, y1=hsrc.male15to24.prev$upperCI, col=1, lty=1, lwd=1.5)
mtext("Men, 15 to 24 y", 3, -1.2, adj=0.05, font=2, cex=0.8)
##
plot(1970:2012, rowMeans(male25to49prev.rmat1), type="l", lty=1, lwd=1.5, col=4, xlim=c(1995, 2012), ylim=c(0, 0.25), ylab="", xlab="")
matlines(1970:2012, t(apply(male25to49prev.rmat1, 1, quantile, c(0.025, 0.975))), lty=2, lwd=1.5, col=4)
points(sa.ageprev$years, hsrc.male25to49.prev$mean, lwd=1.5, cex=0.8)
segments(sa.ageprev$years, y0=hsrc.male25to49.prev$lowerCI, y1=hsrc.male25to49.prev$upperCI, col=1, lty=1, lwd=1.5)
mtext("Men, 25 to 49 y", 3, -1.2, adj=0.05, font=2, cex=0.8)
##
plot(1970:2012, rowMeans(male50plusprev.rmat1), type="l", lty=1, lwd=1.5, col=4, xlim=c(1995, 2012), ylim=c(0, 0.25), ylab="", xlab="")
matlines(1970:2012, t(apply(male50plusprev.rmat1, 1, quantile, c(0.025, 0.975))), lty=2, lwd=1.5, col=4)
points(sa.ageprev$years, hsrc.male50plus.prev$mean, lwd=1.5, cex=0.8)
segments(sa.ageprev$years, y0=hsrc.male50plus.prev$lowerCI, y1=hsrc.male50plus.prev$upperCI, col=1, lty=1, lwd=1.5)
mtext("Men, 50+ y", 3, -1.2, adj=0.05, font=2, cex=0.8)
##
plot(1970:2012, rowMeans(female15to24prev.rmat1), type="l", lty=1, lwd=1.5, col=4, xlim=c(1995, 2012), ylim=c(0, 0.25), ylab="", xlab="")
matlines(1970:2012, t(apply(female15to24prev.rmat1, 1, quantile, c(0.025, 0.975))), lty=2, lwd=1.5, col=4)
points(sa.ageprev$years, hsrc.female15to24.prev$mean, lwd=1.5, cex=0.8)
segments(sa.ageprev$years, y0=hsrc.female15to24.prev$lowerCI, y1=hsrc.female15to24.prev$upperCI, col=1, lty=1, lwd=1.5)
mtext("Women, 15 to 24 y", 3, -1.2, adj=0.05, font=2, cex=0.8)
##
plot(1970:2012, rowMeans(female25to49prev.rmat1), type="l", lty=1, lwd=1.5, col=4, xlim=c(1995, 2012), ylim=c(0, 0.4), ylab="", xlab="")
matlines(1970:2012, t(apply(female25to49prev.rmat1, 1, quantile, c(0.025, 0.975))), lty=2, lwd=1.5, col=4)
points(sa.ageprev$years, hsrc.female25to49.prev$mean, lwd=1.5, cex=0.8)
segments(sa.ageprev$years, y0=hsrc.female25to49.prev$lowerCI, y1=hsrc.female25to49.prev$upperCI, col=1, lty=1, lwd=1.5)
mtext("Women, 25 to 49 y", 3, -1.2, adj=0.05, font=2, cex=0.8)
##
plot(1970:2012, rowMeans(female50plusprev.rmat1), type="l", lty=1, lwd=1.5, col=4, xlim=c(1995, 2012), ylim=c(0, 0.15), ylab="", xlab="")
matlines(1970:2012, t(apply(female50plusprev.rmat1, 1, quantile, c(0.025, 0.975))), lty=2, lwd=1.5, col=4)
points(sa.ageprev$years, hsrc.female50plus.prev$mean, lwd=1.5, cex=0.8)
segments(sa.ageprev$years, y0=hsrc.female50plus.prev$lowerCI, y1=hsrc.female50plus.prev$upperCI, col=1, lty=1, lwd=1.5)
mtext("Women, 50+ y", 3, -1.2, adj=0.05, font=2, cex=0.8)
##


quartz(h=4.5, w=8, pointsize=10)

par(mfrow=c(2, 3), mar=c(2, 2.5, 0.5, 0.5), mgp=c(2, 0.5, 0), tcl=-0.25, cex=1, las=1)
##
plot(1970:2012, rowMeans(male15to24prev.incrr1), type="l", lty=1, lwd=1.5, col=2, xlim=c(1995, 2012), ylim=c(0, 0.15), ylab="", xlab="")
matlines(1970:2012, t(apply(male15to24prev.incrr1, 1, quantile, c(0.025, 0.975))), lty=2, lwd=1.5, col=2)
points(sa.ageprev$years, hsrc.male15to24.prev$mean, lwd=1.5, cex=0.8)
segments(sa.ageprev$years, y0=hsrc.male15to24.prev$lowerCI, y1=hsrc.male15to24.prev$upperCI, col=1, lty=1, lwd=1.5)
mtext("Men, 15 to 24 y", 3, -1.2, adj=0.05, font=2, cex=0.8)
##
plot(1970:2012, rowMeans(male25to49prev.incrr1), type="l", lty=1, lwd=1.5, col=2, xlim=c(1995, 2012), ylim=c(0, 0.25), ylab="", xlab="")
matlines(1970:2012, t(apply(male25to49prev.incrr1, 1, quantile, c(0.025, 0.975))), lty=2, lwd=1.5, col=2)
points(sa.ageprev$years, hsrc.male25to49.prev$mean, lwd=1.5, cex=0.8)
segments(sa.ageprev$years, y0=hsrc.male25to49.prev$lowerCI, y1=hsrc.male25to49.prev$upperCI, col=1, lty=1, lwd=1.5)
mtext("Men, 25 to 49 y", 3, -1.2, adj=0.05, font=2, cex=0.8)
##
plot(1970:2012, rowMeans(male50plusprev.incrr1), type="l", lty=1, lwd=1.5, col=2, xlim=c(1995, 2012), ylim=c(0, 0.25), ylab="", xlab="")
matlines(1970:2012, t(apply(male50plusprev.incrr1, 1, quantile, c(0.025, 0.975))), lty=2, lwd=1.5, col=2)
points(sa.ageprev$years, hsrc.male50plus.prev$mean, lwd=1.5, cex=0.8)
segments(sa.ageprev$years, y0=hsrc.male50plus.prev$lowerCI, y1=hsrc.male50plus.prev$upperCI, col=1, lty=1, lwd=1.5)
mtext("Men, 50+ y", 3, -1.2, adj=0.05, font=2, cex=0.8)
##
plot(1970:2012, rowMeans(female15to24prev.incrr1), type="l", lty=1, lwd=1.5, col=2, xlim=c(1995, 2012), ylim=c(0, 0.25), ylab="", xlab="")
matlines(1970:2012, t(apply(female15to24prev.incrr1, 1, quantile, c(0.025, 0.975))), lty=2, lwd=1.5, col=2)
points(sa.ageprev$years, hsrc.female15to24.prev$mean, lwd=1.5, cex=0.8)
segments(sa.ageprev$years, y0=hsrc.female15to24.prev$lowerCI, y1=hsrc.female15to24.prev$upperCI, col=1, lty=1, lwd=1.5)
mtext("Women, 15 to 24 y", 3, -1.2, adj=0.05, font=2, cex=0.8)
##
plot(1970:2012, rowMeans(female25to49prev.incrr1), type="l", lty=1, lwd=1.5, col=2, xlim=c(1995, 2012), ylim=c(0, 0.4), ylab="", xlab="")
matlines(1970:2012, t(apply(female25to49prev.incrr1, 1, quantile, c(0.025, 0.975))), lty=2, lwd=1.5, col=2)
points(sa.ageprev$years, hsrc.female25to49.prev$mean, lwd=1.5, cex=0.8)
segments(sa.ageprev$years, y0=hsrc.female25to49.prev$lowerCI, y1=hsrc.female25to49.prev$upperCI, col=1, lty=1, lwd=1.5)
mtext("Women, 25 to 49 y", 3, -1.2, adj=0.05, font=2, cex=0.8)
##
plot(1970:2012, rowMeans(female50plusprev.incrr1), type="l", lty=1, lwd=1.5, col=2, xlim=c(1995, 2012), ylim=c(0, 0.15), ylab="", xlab="")
matlines(1970:2012, t(apply(female50plusprev.incrr1, 1, quantile, c(0.025, 0.975))), lty=2, lwd=1.5, col=2)
points(sa.ageprev$years, hsrc.female50plus.prev$mean, lwd=1.5, cex=0.8)
segments(sa.ageprev$years, y0=hsrc.female50plus.prev$lowerCI, y1=hsrc.female50plus.prev$upperCI, col=1, lty=1, lwd=1.5)
mtext("Women, 50+ y", 3, -1.2, adj=0.05, font=2, cex=0.8)
##
