library(mvtnorm)

logit <- function(p){ return(log(p / (1-p))) }
invlogit <- function(x){ return(exp(x) / (1+exp(x))) }

create.Rmat <- function(m.mean, m.sd, f.mean, f.sd, corr, margdist = "gamma", copula = "gaussian", VERSION="C"){

  if(VERSION != "R"){
    if(margdist != "gamma" || copula != "gaussian")
      stop("not implemented")

    return(.Call("createRmatR", m.mean, m.sd, f.mean, f.sd, corr))
  }

  ##################################################################################

  
  Rmat <- matrix(0, AG, AG)
  
  ## get the quantiles of the marginal distribution
  if(margdist == "gamma"){
    m.rate <- m.mean / m.sd^2
    m.shape <- m.mean^2 / m.sd^2

    f.rate <- f.mean / f.sd^2
    f.shape <- f.mean^2 / f.sd^2

    m.quant <- qnorm(pgamma(c(0:(AG-4)*5, Inf), m.shape, m.rate))
    f.quant <- qnorm(pgamma(c(0:(AG-4)*5, Inf), f.shape, f.rate))
  }


  ## generate discretized PMF based on copula
  if(copula == "gaussian"){
    Sig <- matrix(c(1, corr, corr, 1), 2)
    Rmat[age15plus.idx, age15plus.idx] <- outer(1:(AG-3), 1:(AG-3),
                                                Vectorize(function(i, j) pmvnorm(c(m.quant[i], f.quant[j]), c(m.quant[i+1], f.quant[j+1]), sigma=Sig)))
  }

  return(Rmat)
}

create.incrr <- function(m.mean, m.sd, f.mean, f.sd, f.incrr, param, margdist = "gamma"){

  param$inc.agerat.male <- numeric(AG)
  param$inc.agerat.female <- numeric(AG)
  param$inc.sexrat <- f.incrr
  
  if(margdist == "gamma"){
    m.rate <- m.mean / m.sd^2
    m.shape <- m.mean^2 / m.sd^2

    f.rate <- f.mean / f.sd^2
    f.shape <- f.mean^2 / f.sd^2

    m.cum <- pgamma(c(0:(AG-4)*5, Inf), m.shape, m.rate)
    f.cum <- pgamma(c(0:(AG-4)*5, Inf), f.shape, f.rate)

    param$inc.agerat.male[4:AG] <- m.cum[-1] - m.cum[-(AG-2)]
    param$inc.agerat.female[4:AG] <- f.cum[-1] - f.cum[-(AG-2)]
  }

  return(param)  # function returns updated param object
}


## age/sex-specific HIV prevalence incorporating transmission reduction for those on ART
fnEffectiveAgePrev <- function(X) {(rowSums(X[,,-1,1],,2) + relinfect.ART*rowSums(X[,,-1,-1],,2))/rowSums(X,,2)}  


fnAgeInc <- function(X, r, iota, param, year.idx){

  if(!is.null(param$inc.agerat.male)){
    inc.rate.15to49 <-  r * ((sum(X[,a15to49.idx,-1,1]) + relinfect.ART*sum(X[,a15to49.idx,-1,-1]))/sum(X[,a15to49.idx,,]) + iota)
    inc.rr <- rbind(param$inc.agerat.male, param$inc.sexrat*param$inc.agerat.female)
    age.inc <- inc.rr*inc.rate.15to49/(sum(X[,a15to49.idx,1,1]*inc.rr[,a15to49.idx])/sum(X[,a15to49.idx,1,1]))

  } else if(!is.null(param$Rmat)){
    eff.ageprev <- fnEffectiveAgePrev(X)
    age.inc <- rbind(r * (eff.ageprev[2,] + iota) %*% t(param$Rmat), # incidence rate in men, function of prevalence in women
                     r * (eff.ageprev[1,] + iota) %*% param$Rmat)    # incidence rate in women, function of prevalence in men
  } else {
    stop(paste("'", param$model, "' model is not implemented.", sep=""))
  }
  
  return(age.inc)
}
