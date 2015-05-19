
dir <- getwd()
## setwd("~/anclik/")
setwd("~/Documents/Code/R/anclik/")
source("anclik.R")
setwd(dir)
rm(dir)

source("R/epp.R")
source("R/spectrum.R")
source("R/generics.R")

source("R/IMIS.R")


#################
####  Prior  ####
#################

logiota.unif.prior <- c(log(1e-14), log(0.0025))
tau2.prior.rate <- 0.5

invGammaParameter <- 0.001   #Inverse gamma parameter for tau^2 prior for spline
ancbias.pr.mean <- 0.15
ancbias.pr.sd <- 1.0
muSS <- 1/11.5               #1/duration for r steady state prior

lprior <- function(theta, fp){

  nk <- fp$numKnots
  tau2 <- exp(theta[nk+3])

  return(sum(dnorm(theta[3:nk], 0, sqrt(tau2), log=TRUE)) +
         dunif(theta[nk+1], logiota.unif.prior[1], logiota.unif.prior[2], log=TRUE) + 
         dnorm(theta[nk+2], ancbias.pr.mean, ancbias.pr.sd, log=TRUE) +
         log(densigamma(tau2, invGammaParameter, invGammaParameter)))         
}

################################
####                        ####
####  HH survey likelihood  ####
####                        ####
################################

library(mvtnorm)
library(pscl)  # for densigamma()


fnHHSll <- function(qM, hhslik.dat){
  return(sum(dnorm(hhslik.dat$W.hhs, qM[hhslik.dat$idx], hhslik.dat$sd.W.hhs, log=TRUE)))
}


fnCreateParam <- function(theta, fp){

  u <- theta[1:fp$numKnots]
  beta <- numeric(fp$numKnots)
  beta[1] <- u[1]
  beta[2] <- u[1]+u[2]
  for(i in 3:fp$numKnots)
    beta[i] <- -beta[i-2] + 2*beta[i-1] + u[i]
    
  return(list(rvec = as.vector(fp$rvec.spldes %*% beta),
              iota = exp(theta[fp$numKnots+1]),
              ancbias = theta[fp$numKnots+2]))
}

ll <- function(theta, fp, likdat){
  theta.last <<- theta
  param <- fnCreateParam(theta, fp)

  if(min(param$rvec)<0 || max(param$rvec)>20) # Test positivity of rvec
    return(-Inf) 
  
  if(inherits(fp, "specfp"))  ## TODO: revise to use generic functions
    mod <- fnSpectrum(param, fp)
  else
    mod <- fnEPP(param, fp)

  qM.all <- qnorm(prev(mod))
  qM.preg <- if(fp$pregprev) qnorm(fnPregPrev(mod, fp)) else qM.all

  if(any(is.na(qM.all)) || any(qM.all[-1] == -Inf))
    return(-Inf)

  ll.anc <- log(fnANClik(qM.preg+param$ancbias, likdat$anclik.dat))
  ll.hhs <- fnHHSll(qM.all, likdat$hhslik.dat)

  if(exists("equil.rprior", where=fp) && fp$equil.rprior){
    rvec.ann <- param$rvec[fp$proj.steps %% 1 == 0.5]
    equil.rprior.mean <- muSS/(1-pnorm(qM.all[likdat$lastdata.idx]))
    equil.rprior.sd <- sqrt(mean((muSS/(1-pnorm(qM.all[lastdata.idx - 10:1])) - rvec.ann[lastdata.idx - 10:1])^2))  # empirical sd based on 10 previous years
    ll.rprior <- sum(dnorm(rvec.ann[lastdata.idx:length(qM.all)], equil.rprior.mean, equil.rprior.sd, log=TRUE))  # prior starts in last data year (Dan's starts in next year)
  } else
    ll.rprior <- 0
  
  return(ll.anc+ll.hhs+ll.rprior)
}


##########################
####  IMIS functions  ####
##########################

## Note: requires fp and likdat to be in global environment

sample.prior <- function(n){
  
  mat <- matrix(NA, n, fp$numKnots+3)
  
  ## sample penalty variance
  tau2 <- rexp(n, tau2.prior.rate)                  # variance of second-order spline differences
  
  mat[,1] <- rnorm(n, 1.5, 1)                                                     # u[1]
  mat[,2:fp$numKnots] <- rnorm(n*(fp$numKnots-1), 0, sqrt(tau2))                  # u[2:numKnots]
  mat[,fp$numKnots+1] <-  runif(n, logiota.unif.prior[1], logiota.unif.prior[2])  # iota
  mat[,fp$numKnots+2] <-  rnorm(n, ancbias.pr.mean, ancbias.pr.sd)                # ancbias parameter
  mat[,fp$numKnots+3] <- log(tau2)                                                # tau2
  
  
  return(mat)
}

prior <- function(theta){
  return(unlist(lapply(seq_len(nrow(theta)), function(i) return(exp(lprior(theta[i,], fp))))))
}

likelihood <- function(theta){
  return(unlist(lapply(seq_len(nrow(theta)), function(i) return(exp(ll(theta[i,], fp, likdat))))))
}


#########################################################
####  Define fixed parameters for fitting scenarios  ####
#########################################################

epp.fp <- list(pregprev=FALSE)
spec.fp <- list(pregprev=FALSE)
## preg.fp <- list("cd4stage.weight" = c(1.3, 0.6, 0.1, 0.1, 0.0, 0.0, 0.0),  ## for EPP stage weighted prev
##                 "art1yr.weight" = 0.3, pregprev=TRUE)
asfr.fp <- list("age.fertrat" = rep(1, 7), "stage.fertrat" = rep(1, 7), "art.fertrat" = 1, pregprev=TRUE)
age.fp <- list("age.fertrat" = c(1.20, 0.76, 0.71, 0.65, 0.59, 0.53, 0.47), "stage.fertrat" = rep(1, 7), "art.fertrat" = 1, pregprev=TRUE)
stage.fp <- list("age.fertrat" = c(1.59, 1.06, 1.06, 1.00, 0.98, 0.91, 0.86),
                 "stage.fertrat" = c(1.25, 0.6, 0.4, 0.4, 0.3, 0.3, 0.3),
                 "art.fertrat" = 0.6, pregprev=TRUE)
