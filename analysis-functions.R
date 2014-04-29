############################
####  Spline functions  ####
############################

proj.steps <- seq(1, 43.5, 0.1)  # 1970:2012
t0 <- 6
dt <- 0.1
numSplineSteps <- sum(proj.steps >= t0)
numSplines <- 7

library(splines)

fnBSpline <- function(u, q = numSplines, m = 2){
  x <- seq(1,numSplineSteps)

  k <- seq(min(x),max(x),length=q-m)  # knots
  dk <- k[2]-k[1]
  k <- c(k[1]-dk*((m+1):1),k,k[q-m]+dk*(1:(m+1)))

  X <- splineDesign(k,x,ord=m+2)      # model matrix


  b <- u                              # Now the coefficients
  b[1] <- u[1]
  b[2] <- u[2]
  for(i in 3:length(b)){
    b[i] <- 2*b[i-1] - b[i-2] + u[i]
  }

  ## Predicted r
  r <- c(rep(0, sum(proj.steps < t0)), X%*%b)

  return(r)
}


##############################
####  Analysis functions  ####
##############################

fert.idx <- 4:10
a15to49.idx <- 4:10
a15to24.idx <- 4:5
a15plus.idx <- 4:AG

m.idx <- 1
f.idx <- 2

fnRVec <- function(theta){return(fnBSpline(theta[-(1:2)]))}
fnIota <- function(theta){return(exp(theta[1]))}

fnMod <- function(theta){return(fnSpectrum(fnIota(theta), fnRVec(theta)))}

fnASFRadjPrev <- function(mod){
  prev.by.age <- rowSums(mod[,2, fert.idx, -1,],,2)/rowSums(mod[,2, fert.idx,,],,2)
  births.by.age <- t(asfr[,1:dim(mod)[1]]) * rowSums(mod[,2, fert.idx,,],,2)
  return(rowSums(prev.by.age * births.by.age) / rowSums(births.by.age))
}

fnANCprev <- function(mod){
  hivn.by.age <- rowSums(mod[,2, fert.idx, 1,],,2)
  hivp.by.age <- rowSums(mod[,2, fert.idx, -1,],,2)
  births.by.age <- t(asfr[,1:dim(mod)[1]]) * (hivn.by.age + hivp.by.age)
  frac.hivp.moth <- 1.0 - hivn.by.age/(sweep(hivp.by.age, 2, fert.rat, "*")+hivn.by.age)
  return(rowSums(births.by.age * frac.hivp.moth) / rowSums(births.by.age))
}

prev <- function(mod, age.idx=a15to49.idx, sex.idx=c(m.idx, f.idx)) return(rowSums(mod[,sex.idx, age.idx,-1,])/rowSums(mod[,sex.idx, age.idx,,]))

## fn15to49inc <- function(mod)

## inc.rate.15to49 <- rVec[ts == proj.steps] * ((sum(X[,age15to49.idx,-1,1]) + relinfect.ART*sum(X[,age15to49.idx,-1,-1]))/sum(X[,age15to49.idx,,]) + (ts == t0)*iota)
