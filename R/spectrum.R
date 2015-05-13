####################################
####                            ####
####  Load the C share library  ####
####                            ####
####################################

dyn.load("C++/rlib.so")


#############################
####                     ####
####  Declare constants  ####
####                     ####
#############################

## TODO: Put this in a namespace

"incr<-" <- function(x, value) { x + value } # increment operator, from Hmisc

NG <- 2   # number of sexes
AG <- 17  # number of age groups
DS <- 8   # number of disease stages
TS <- 4   # ART treatment duration stages

m.idx <- 1
f.idx <- 2

fert.idx <- 4:10
age15to49.idx <- 4:10
age15plus.idx <- 4:AG


#########################################
####                                 ####
####  Create fixed model parameters  ####
####                                 ####
#########################################

library(splines)

CreateSpectrumFixpar <- function(projp, demp, dt = 0.1, proj.start = projp$yr.start+dt*ceiling(1/(2*dt)),
                                 proj.end = projp$yr.end+dt*ceiling(1/(2*dt)),
                                 relinfectART = 0.3, time.epi.start = proj.start) {

  proj.steps <- seq(proj.start, proj.end, dt)

  ## Function to linearly interpolate rates that are annual inputs to each timestep
  interpolate.ts <- function(arr, xout=proj.steps, dm=1){
    apply(arr, seq_len(length(dim(arr)))[-dm],
          function(y) approx(x=as.numeric(dimnames(arr)[[dm]])+0.5,  # +0.5 for mid-year
                             y=y, xout=xout)$y)
  }

  ## Parameters that change over time will have timestep as the last index
  ## (e.g. [SEX, AGE, TS]), so that values for each ts will be in a contigous
  ## block of memory. aperm() does this.

  init.pop <- interpolate.ts(demp$basepop, proj.start)  # linearly interpolate between basepop years

  mx.ts <- aperm(interpolate.ts(demp$mx), c(2:3,1))
  srb.ts <- approx(as.numeric(names(demp$srb)), demp$srb/100, proj.steps)$y
  asfr.ts <- aperm(interpolate.ts(demp$asfr), 2:1)

  agesex.incrr <- projp$inc.agerat
  agesex.incrr[, f.idx,] <- sweep(agesex.incrr[, f.idx,], 1, projp$inc.sexrat, "*")
  agesex.incrr.ts <- aperm(interpolate.ts(agesex.incrr), c(2:3,1))


  artelig.idx <- match(projp$arteligthresh.15plus, c(1, 2, 500, 350, 250, 200, 100, 50))
  artelig.idx.ts <- approx(as.numeric(names(projp$arteligthresh.15plus)), artelig.idx, proj.steps, "constant", rule=2)$y


  ## Number of persons who should be on ART at end of timestep
  ## Offset by 1 year because Spectrum number on ART are Dec 31
  ## !! Currently only dealing with numbers on ART, assuming no increase in coverage after last number
  ## !! Doesn't handle percentage coverage inputs
  m.artnum.15plus <- projp$artnum.15plus[!projp$artnumperc.15plus[,m.idx], m.idx]
  f.artnum.15plus <- projp$artnum.15plus[!projp$artnumperc.15plus[,f.idx], f.idx]
  artnum.15plus.ts <- rbind("Male"  =approx(as.numeric(names(m.artnum.15plus))+1-dt, m.artnum.15plus, proj.steps, rule=2)$y,
                            "Female"=approx(as.numeric(names(f.artnum.15plus))+1-dt, f.artnum.15plus, proj.steps, rule=2)$y)


  ###########################
  ##  r-spline parameters  ##
  ###########################

  numKnots <- 7
  proj.dur <- diff(range(proj.steps))
  rvec.knots <- seq(min(proj.steps) - 3*proj.dur/(numKnots-3), max(proj.steps) + 3*proj.dur/(numKnots-3), proj.dur/(numKnots-3))
  rvec.spldes <- splineDesign(rvec.knots, proj.steps)


  fp <- list(dt               = dt,
             proj.steps       = proj.steps,
             ts.epi.start     = as.integer(which(proj.steps == time.epi.start)),
             init.pop         = init.pop,
             mx.ts            = mx.ts,
             asfr.ts          = asfr.ts,
             srb.ts           = srb.ts,
             agesex.incrr.ts  = agesex.incrr.ts,
             age.fertrat      = projp$fert.rat,
             stage.fertrat    = c(1.25, 0.6, 0.4, 0.4, 0.3, 0.3, 0.3),
             art.fertrat      = 0.6,
             vert.trans       = 0,
             cd4.prog         = projp$cd4.prog,
             cd4.initdist     = projp$cd4.initdist,
             cd4.art.mort     = projp$cd4.art.mort,
             relinfectART     = relinfectART,
             artelig.idx.ts   = as.integer(artelig.idx.ts),
             artnum.15plus.ts = artnum.15plus.ts,
             numKnots         = numKnots,
             rvec.spldes      = rvec.spldes)
  class(fp) <- "specfp"

  return(fp)
}

update.specfp <- function(fp, ..., keep.attr=TRUE){
  ## adapted from 'survey' package: https://github.com/cran/survey/blob/master/R/survey.R

  dots<-substitute(list(...))[-1]
  newnames<-names(dots)
  
  for(j in seq(along=dots)){
    if(keep.attr)
      attr <- attributes(fp[[newnames[j]]])
    fp[[newnames[j]]]<-eval(dots[[j]],fp, parent.frame())
    if(keep.attr)
      attributes(fp[[newnames[j]]]) <- attr
  }
  
  return(fp)
}


##########################
####                  ####
####  Spectrum model  ####
####                  ####
##########################

fnSpectrum <- function(param, fp, VERSION = "C"){

  
  if(VERSION != "R"){
    OUT.STEPS <- sum(fp$proj.steps %% 1 == fp$dt*floor((1/fp$dt)/2))
    mod <- .Call("fnSpectrumR", param, fp, OUT.STEPS)
    class(mod) <- "spec"
    return(mod)
  }

  
  ##################################################################################

  dt <- fp$dt
  proj.steps <- fp$proj.steps

  ## initialize output
  X.out <- array(NA, c(sum(proj.steps %% 1 == dt*floor((1/dt)/2)), NG, AG, DS, TS))  # TODO: permute with outyear last dim

  ## initialize population
  X <- array(0, c(NG, AG, DS, TS))
  X[,,1,1] <- fp$init.pop

  for(ts in 1:length(proj.steps)){

    if(proj.steps[ts] %% 1 == dt*floor((1/dt)/2)) # store the result mid-year (or ts before mid-year)
      X.out[ceiling(ts*dt),,,,] <- X

    grad <- array(0, c(NG, AG, DS, TS))

    ## ageing
    incr(grad[,-AG,,]) <- -X[,-AG,,]*0.2  # ageing
    incr(grad[,-1,,]) <- X[,-AG,,]*0.2    # ageing

    ## natural mortality
    incr(grad) <- -sweep(X, 1:2, fp$mx.ts[,,ts], "*")

    ## fertility
    births.by.age <- rowSums(X[f.idx,fert.idx,,])*fp$asfr.ts[,ts]
    hivn.by.age <- rowSums(X[f.idx, fert.idx, 1,])
    weighted.hivp.by.age <- (rowSums(sweep(rowSums(X[f.idx, fert.idx, -1,1:3],,2), 2, fp$stage.fertrat, "*")) + fp$art.fertrat * rowSums(X[f.idx, fert.idx,-1,4])) * fp$age.fertrat
    frac.hivp.births <- fp$vert.trans * weighted.hivp.by.age / (weighted.hivp.by.age + hivn.by.age)  # fraction of births HIV+ by age group

    incr(grad[,1,-1,1]) <- sum(births.by.age*frac.hivp.births) * c(fp$srb[ts], 1.0)/(fp$srb[ts]+1.0) * fp$cd4.initdist[,1,] # HIV+ children born
    incr(grad[,1,1,1]) <- sum(births.by.age*(1.0 - frac.hivp.births)) * c(fp$srb.ts[ts], 1.0)/(fp$srb.ts[ts]+1.0) # HIV- children born


    if (ts >= fp$ts.epi.start){

      ## incidence
      ## age.inc <- fnAgeInc(X, param$rVec[ts == proj.steps], (ts == t0)*param$iota, param, year.idx)
      incrate.15to49 <- param$rvec[ts] * (sum(X[,age15to49.idx,-1,1]) + fp$relinfectART*sum(X[,age15to49.idx,-1,-1]))/sum(X[,age15to49.idx,,]) + param$iota*(ts == fp$ts.epi.start)
      inc.rr <- fp$agesex.incrr.ts[,,ts]
      age.inc <- inc.rr*incrate.15to49/(sum(X[,age15to49.idx,1,1]*inc.rr[,age15to49.idx])/sum(X[,age15to49.idx,1,1]))

      incr(grad[,,1,1]) <- -age.inc*X[,,1,1]  # remove incident infections from susceptibles
      incr(grad[,,-1,1]) <- sweep(fp$cd4.initdist, 1:2, age.inc*X[,,1,1], "*")  # add incident infections

      ## print(paste("ts", ts-1, "Xhivn", sum(X[,age15to49.idx,1,1]), "Xhivp_noart", sum(X[,age15to49.idx,-1,1]), "Xart", sum(X[,age15to49.idx,,-1]), "incrate", incrate.15to49))

      ## disease progression and mortality
      incr(grad[,,-c(1,DS),1]) <- -fp$cd4.prog * X[,,-c(1, DS),1]  # remove cd4 stage progression (untreated)
      incr(grad[,,-(1:2),1]) <- fp$cd4.prog * X[,,-c(1, DS),1]     # add cd4 stage progression (untreated)
      incr(grad[,,-1, 2:3]) <- -2 * X[,,-1, 2:3]                   # remove ART duration progression (HARD CODED 6 months duration)
      incr(grad[,,-1, 3:4]) <- 2 * X[,,-1, 2:3]                    # add ART duration progression (HARD CODED 6 months duration)
      incr(grad[,,-1,]) <- -fp$cd4.art.mort * X[,,-1,]             # HIV mortality

      ## ART initiation
      if (sum(fp$artnum.15plus.ts[,ts]) > 0){

        ts.elig.idx <- fp$artelig.idx.ts[ts]:DS
        artnum.15plus.curr <- rowSums(X[,age15plus.idx,,-1])
        art.15plus.anninits <- (fp$artnum.15plus.ts[,ts] - artnum.15plus.curr)/dt - rowSums(grad[,age15plus.idx,,-1]) # desired change rate minus current exits
        artelig.15plus <- X[, age15plus.idx, ts.elig.idx, 1]
        expect.mort.weight <- sweep(fp$cd4.art.mort[, age15plus.idx, ts.elig.idx - 1, 1], 1,
                                    rowSums(artelig.15plus * fp$cd4.art.mort[, age15plus.idx, ts.elig.idx - 1, 1]), "/")
        artinit.weight <- sweep(expect.mort.weight, 1, 1/rowSums(artelig.15plus), "+")/2
        artinit.ann.ts <- pmin(sweep(artinit.weight * artelig.15plus, 1, art.15plus.anninits, "*"),
                               artelig.15plus/dt,
                               artelig.15plus/dt + grad[, age15plus.idx, ts.elig.idx,1])

        incr(grad[,age15plus.idx, ts.elig.idx, 1]) <- -artinit.ann.ts
        incr(grad[,age15plus.idx, ts.elig.idx, 2]) <- artinit.ann.ts

      }  # if (sum(fp$artnum.15plus.ts[,ts]) > 0)
    }  # if (ts >= fp$ts.epi.start)

    ## do projection (euler integration) ##
    incr(X) <- dt*grad
  } # for(ts in 1:length(proj.steps))

  class(X.out) <- "spec"
  return(X.out)
}


##############################
####  Analysis functions  ####
##############################

prev.spec <- function(mod, age.idx=age15to49.idx, sex.idx=c(m.idx, f.idx)){
  return(rowSums(mod[,sex.idx, age.idx,-1,])/rowSums(mod[,sex.idx, age.idx,,]))
}

fnPregPrev.spec <- function(mod, fp){
  births.by.age <- t(fp$asfr.ts[,fp$proj.steps %% 1 == 0.5]) * rowSums(mod[,f.idx, fert.idx,,],,2)
  hivn.by.age <- rowSums(mod[,f.idx, fert.idx, 1,],,2)
  weighted.hivp.by.age <- sweep(rowSums(sweep(rowSums(mod[,f.idx, fert.idx, -1,1:3],,3), 3, fp$stage.fertrat, "*"),,2) + fp$art.fertrat * rowSums(mod[,f.idx, fert.idx,-1,4],,2), 2, fp$age.fertrat, "*")
  frac.hivp.moth <- weighted.hivp.by.age/(weighted.hivp.by.age+hivn.by.age)
  return(rowSums(births.by.age * frac.hivp.moth) / rowSums(births.by.age))
}

incid.spec <- function(mod, param, fp, age.idx=age15to49.idx, sex.idx=c(m.idx, f.idx)){
  param$rvec[fp$proj.steps %% 1 == 0.5] * (rowSums(mod[, , age15to49.idx, -1, 1]) + fp$relinfectART * rowSums(mod[, , age15to49.idx, -1, -1])) / rowSums(mod[,, age15to49.idx,1,])
}
