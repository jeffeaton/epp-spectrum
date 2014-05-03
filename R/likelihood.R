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

source("C++/functions.R")
source("analysis-functions.R")


################
####  Data  ####
################

source("../south-africa/sa-prevalence-data.R")
sa.anc$idx <- 21:42
sa.15to49$idx <- c(33, 36, 39, 43)

anc.dat <- sa.anc
hhsurv.dat <- sa.15to49


######################
####  Likelihood  ####
######################

library(multicore)
library(pscl)

## Prior parameters
logiota.unif.prior <- c(log(1e-13), log(0.0025))
# tau2.prior.rate <- 0.5
tau2.prior.rate <- 1.0
invGammaParameter <- 0.001   #Inverse gamma parameter for tau^2 prior for spline
## muSS <- 1/11.5               #1/duration for r steady state prior


lprior <- function(theta){
  tau2 <- exp(theta[2])
  return(dunif(theta[1], logiota.unif.prior[1], logiota.unif.prior[2], log=TRUE) +
           log(densigamma(tau2, invGammaParameter, invGammaParameter)) +
           sum(dnorm(theta[5:length(theta)], 0, sqrt(tau2), log=TRUE)))
}

prior <- function(theta){
  if(is.vector(theta))
    theta <- matrix(theta, 1)

  ## return(exp(apply(theta, 1, lprior)))
  nr <- nrow(theta)
  return(unlist(mclapply(1:nr, function(i) return(exp(lprior(theta[i,]))))))
}

sample.prior <- function(n){

 mat <- matrix(NA, n, numSplines+2)
 mat[,1] <- runif(n, logiota.unif.prior[1], logiota.unif.prior[2])
 tau2 <- rexp(n, tau2.prior.rate)  ## !! NOT SURE WHY THIS IS SAMPLING FROM AN EXPONENTIAL
 mat[,2] <- log(tau2)
 mat[,3] <- rnorm(n, 1, 2) ## NOT SURE WHAT THIS IS
 mat[,-(1:3)] <- rnorm(n*(numSplines-1), 0, sqrt(tau2))

 return(mat)
}

ll <- function(theta){

  iota <- exp(theta[1])  # initial epidemic pulse
  rVec <- fnBSpline(theta[-(1:2)])

  # if(min(rVec)<0 || max(rVec)>20) # Test positivity constraint
  if(min(rVec)<0 || max(rVec)>3) # Test positivity constraint
    return(-Inf) 
  
  mod <- fnSpectrum(iota, rVec)
  lM15to49 <- logit(prev(mod)[hhsurv.dat$idx])
  if(any(is.na(lM15to49)) || any(lM15to49 > logit(0.95)))
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
  ll.hhsurv <- sum(-log(hhsurv.dat$logit.var)/2 - (hhsurv.dat$logit.p - lM15to49)^2/(2*hhsurv.dat$logit.var))

  if(any(is.na(ll.anc)) || any(is.na(ll.hhsurv)))
    return(-Inf)

  return(ll.anc + ll.hhsurv)
}

likelihood <- function(theta){
  if(is.vector(theta))
     theta <- matrix(theta, 1)
  nr <- nrow(theta)
  return(unlist(mclapply(1:nr, function(i) return(exp(ll(theta[i,]))))))
}
