library(gdata)
fert.dist <- read.xls("../south-africa/south-africa-demproj.xlsx", "ASFR", row.names=1, as.is=TRUE)
fert.dist <- data.frame(lapply(fert.dist, as.numeric))
tfr <- as.numeric(read.xls("../south-africa/south-africa-demproj.xlsx", "TFR", row.names=1))
asfr <- mapply("*", fert.dist[1:7,]/500, as.list(tfr))
fert.rat <- read.xls("../south-africa/south-africa-nathist.xlsx", "fertility ratio", row.names=1)$Ratio
rm(fert.dist, tfr)

#################
####  Model  ####
#################

## source("C++/functions.R")
## source("R/spectrum.R")
source("analysis-functions.R")


################
####  Data  ####
################

source("../south-africa/sa-prevalence-data.R")
sa.anc$idx <- 21:42
## sa.15to49$idx <- c(33, 36, 39, 43)
sa.ageprev$idx <- c(33, 36, 39, 43)

anc.dat <- sa.anc
## hhsurv.dat <- sa.15to49
hhsurv.ageprev <- sa.ageprev


######################
####  Likelihood  ####
######################

if(!exists("MODEL")) MODEL <- "rmat"

library(parallel)
options("mc.cores" = 8L)
library(pscl)

## Prior parameters
# logiota.unif.prior <- c(log(1e-14), log(0.0025))
logiota.unif.prior <- c(log(1e-14), log(1e-10))
# tau2.prior.rate <- 0.5
# tau2.prior.rate <- 1.0
tau2.prior.rate <- 0.1
invGammaParameter <- 0.001   #Inverse gamma parameter for tau^2 prior for spline
## muSS <- 1/11.5               #1/duration for r steady state prior


lprior <- function(theta){
  tau2 <- exp(theta[2])
  return(dunif(theta[1], logiota.unif.prior[1], logiota.unif.prior[2], log=TRUE) +
         log(densigamma(tau2, invGammaParameter, invGammaParameter)) +
         sum(dnorm(theta[4 + 1:(numSplines-2)], 0, sqrt(tau2), log=TRUE)) +
         sum(dbeta((theta[2+numSplines+c(2,4)] - 15) / 50, 1.3, 1.3)) +
         sum(dgamma(theta[2+numSplines+c(2,4)], 2.5, 0.05, log=TRUE)) +      # variance of age distribution

         if(MODEL=="rmat")
         dunif(theta[2+numSplines+5], log=TRUE)                             # correlation between m/f age distribution
         else if(MODEL=="incrr")
           dunif(theta[2+numSplines+5], log(0.9), log(1.8), log=TRUE))
         
}

prior <- function(theta){
  if(is.vector(theta))
    theta <- matrix(theta, 1)

  ## return(exp(apply(theta, 1, lprior)))
  nr <- nrow(theta)
  return(unlist(mclapply(1:nr, function(i) return(exp(lprior(theta[i,]))))))
}

sample.prior <- function(n){

 ## mat <- matrix(NA, n, numSplines+2)
 ## mat[,1] <- runif(n, logiota.unif.prior[1], logiota.unif.prior[2])
 ## tau2 <- rexp(n, tau2.prior.rate)
 ## mat[,2] <- log(tau2)
 ## mat[,3] <- rnorm(n, 1, 2)
 ## mat[,-(1:3)] <- rnorm(n*(numSplines-1), 0, sqrt(tau2))

  mat <- matrix(NA, n, 2+numSplines+5)

  mat[,1] <- runif(n, logiota.unif.prior[1], logiota.unif.prior[2])  # iota

  if(MODEL=="rmat")
    tau2 <- rexp(n, tau2.prior.rate)                                   # sd of spline
  else if(MODEL=="incrr")
    tau2 <- rexp(n, 0.5)                                   # sd of spline
  
  mat[,2] <- log(tau2)

  if(MODEL=="rmat")
    mat[,3] <- rnorm(n, 15, 7)
  else if(MODEL == "incrr")
    mat[,3] <- rnorm(n, 1.5, 1)

  mat[,3+1:(numSplines-1)] <- rnorm(n*(numSplines-1), 0, sqrt(tau2))
  mat[,2+numSplines+c(1,3)] <- 15 + 50*rbeta(2*n, 1.3, 1.3)
  mat[,2+numSplines+c(2,4)] <- rgamma(2*n, 2.5, 0.05)

  if(MODEL=="rmat")
    mat[,2+numSplines+5] <- runif(n)
  else if(MODEL=="incrr")
    mat[,2+numSplines+5] <- runif(n, log(0.9), log(1.8))
    
 return(mat)
}

create.param <- function(theta){

  param <- list(iota = exp(theta[1]),  # initial epidemic pulse
                rVec = fnBSpline(theta[2+1:numSplines]))

  if(MODEL=="rmat"){
    Rmat.par <- theta[2+numSplines+1:5]
    Rmat <- create.rmat(Rmat.par[1], Rmat.par[2], Rmat.par[3], Rmat.par[4], Rmat.par[5])
    param$Rmat <- Rmat
  } else if(MODEL=="incrr"){
    par <- theta[2+numSplines+1:5]
    param <- create.incrr(par[1], par[2], par[3], par[4], exp(par[5]), param)
  }

  return(param)
}

ll <- function(theta, ANCPREV_LIKELIHOOD = TRUE, DEBUG = FALSE){

  if(DEBUG)
    lasttheta <<- theta

  param <- create.param(theta)

  # if(min(rVec)<0 || max(rVec)>20) # Test positivity constraint
  if(min(param$rVec)<0 || max(param$rVec)>100) # Test positivity constraint
    return(-Inf) 
  
  mod <- fnSpectrum(param)
  ## lM15to49 <- logit(prev(mod)[hhsurv.dat$idx])
  lMageprev <- logit(ageprev(mod[hhsurv.ageprev$idx,,,,], age15plus.idx, c(1:9, rep(10, 5))))

  
  ## if(any(is.na(lM15to49)) || any(lM15to49 > logit(0.95)))
  ##    return(-Inf)
  if(any(is.na(lMageprev)) || any(lMageprev > logit(0.95)))
    return(-Inf)
  if(ANCPREV_LIKELIHOOD)
    lManc <- logit(fnANCprev(mod)[anc.dat$idx])
  else
    lManc <- logit(prev(mod)[anc.dat$idx])
  if(any(is.na(lManc)) || any(lManc > logit(0.95)))
    return(-Inf)
  
  gamma.unif.pr <- c(0,1)
  S2 <- 1/sum(1/anc.dat$logit.var)
  dbar <- sum((anc.dat$logit.p - lManc)/anc.dat$logit.var)*S2
  d2bar <- sum((anc.dat$logit.p - lManc)^2/anc.dat$logit.var)*S2

  ll.anc <- log(S2)/2 - sum(log(anc.dat$logit.var))/2 + (dbar^2 - d2bar)/S2 + log(pnorm(gamma.unif.pr[2], dbar, sqrt(S2)) - pnorm(gamma.unif.pr[1], dbar, sqrt(S2)))
  ll.hhsurv <- sum(dnorm(lMageprev, hhsurv.ageprev$logit.p, hhsurv.ageprev$logit.se, log=TRUE))

  if(any(is.na(ll.anc)) || any(is.na(ll.hhsurv)))
    return(-Inf)

  return(ll.anc + ll.hhsurv)
}

likelihood <- function(theta, DEBUG=FALSE){
  if(is.vector(theta))
     theta <- matrix(theta, 1)
  nr <- nrow(theta)

  if(DEBUG)
    return(unlist(lapply(1:nr, function(i) return(exp(ll(theta[i,]))))))
  
  return(unlist(mclapply(1:nr, function(i) return(exp(ll(theta[i,]))))))
}
