"incr<-" <- function(x, value) { x + value } # increment operator, from Hmisc

NG <- 2   # number of sexes
AG <- 17  # number of age groups
DS <- 8   # number of disease stages
TS <- 4   # ART treatment duration stages

fert.idx <- 4:10
age15to49.idx <- 4:10
age15plus.idx <- 4:AG

fnSpectrum <- function(param, dt = 0.1){

  ## initialize output
  X.out <- array(NA, c(max(floor(proj.steps - dt*floor((1/dt)/2))),  NG, AG, DS, TS))

  ## initialize population
  X <- array(0, c(NG, AG, DS, TS))
  X[1,,1,1] <- init.pop[1:AG]
  X[2,,1,1] <- init.pop[AG+1:AG]

  for(ts in proj.steps){
    year.idx <- floor(ts)

    if(ts %% 1 == dt*floor((1/dt)/2)) # store the result mid-year (or ts before mid-year)
      X.out[year.idx,,,,] <- X

    grad <- array(0, c(NG, AG, DS, TS))

    ## ageing
    incr(grad[,-AG,,]) <- -X[,-AG,,]*0.2  # ageing
    incr(grad[,-1,,]) <- X[,-AG,,]*0.2    # ageing

    ## natural mortality
    incr(grad[1,,,]) <- -X[1,,,]*mx.male[,year.idx]    # male mortality
    incr(grad[2,,,]) <- -X[2,,,]*mx.female[,year.idx]  # female mortality

    ## fertility
    births.by.age <- rowSums(X[2,fert.idx,,])*asfr[,year.idx]
    frac.hivp.births <- vert.trans*(1.0 - (rowSums(X[2, fert.idx, 1,]))/(fert.rat*rowSums(X[2, fert.idx, -1,]) + rowSums(X[2, fert.idx, 1,]))) # fraction of births HIV+ by age group

    incr(grad[,1,-1,1]) <- sum(births.by.age*frac.hivp.births) * c(srb[year.idx], 100)/(srb[year.idx]+100) * cd4.initdist[,1,] # HIV+ children born
    incr(grad[,1,1,1]) <- sum(births.by.age*(1.0 - frac.hivp.births)) * c(srb[year.idx], 100)/(srb[year.idx]+100) # HIV- children born


    if(ts >= t0){
      ## incidence
      age.inc <- fnAgeInc(X, param$rVec[ts == proj.steps], (ts == t0)*param$iota, param, year.idx)
      incr(grad[,,1,1]) <- -age.inc*X[,,1,1] # remove incident infections from susceptibles
      incr(grad[,,-1,1]) <- array(apply(cd4.initdist, 3, "*", age.inc * X[,,1,1]), c(NG, AG, DS-1)) # add incident infections
      
      ## disease progression and mortality
      incr(grad[,,-c(1,DS),1]) <- -cd4.prog * X[,,-c(1, DS),1]    # remove cd4 stage progression (untreated)
      incr(grad[,,-(1:2),1]) <- cd4.prog * X[,,-c(1, DS),1]       # add cd4 stage progression (untreated)
      incr(grad[,,-1, 2:3]) <- -2 * X[,,-1, 2:3]                  # remove ART duration progression (HARD CODED 6 months duration)
      incr(grad[,,-1, 3:4]) <- 2 * X[,,-1, 2:3]                   # add ART duration progression (HARD CODED 6 months duration)
      incr(grad[,,-1,]) <- - cd4.art.mort * X[,,-1,]              # HIV mortality
      
      ## ART initiation
      numart.15plus.ts <- sum(c(ts %% 1 + dt, 1-(ts %% 1 + dt)) * artnum.15plus[c(year.idx, max(year.idx-1, 1))])
      if(numart.15plus.ts > 0){
        numart.15plus.curr <- sum(X[,age15plus.idx,,-1])
        art.15plus.anninits <- (numart.15plus.ts - numart.15plus.curr)/dt - sum(grad[,age15plus.idx,,-1]) # desired change rate minus current exits

        artelig.15plus <- X[, age15plus.idx, artelig.idx[year.idx]:DS, 1]
        expect.mort.weight <- cd4.art.mort[, age15plus.idx, artelig.idx[year.idx]:DS - 1, 1] / sum(artelig.15plus * cd4.art.mort[, age15plus.idx, artelig.idx[year.idx]:DS - 1, 1])
        artinit.weight <- (expect.mort.weight + 1/sum(artelig.15plus))/2  # average eligibility and expected mortality
        artinit.ann.ts <- art.15plus.anninits * artinit.weight * artelig.15plus
        artinit.ann.ts[artinit.ann.ts > artelig.15plus/dt] <- artelig.15plus[artinit.ann.ts > artelig.15plus/dt]/dt  # check that don't initiate more than the number eligible
      
        incr(grad[,age15plus.idx, artelig.idx[year.idx]:DS, 1]) <- -artinit.ann.ts
        incr(grad[,age15plus.idx, artelig.idx[year.idx]:DS, 2]) <- artinit.ann.ts
      } # if(numart.15plus.ts > 0)
    } # if(ts >= t0)

    ## do projection (euler integration) ##
    incr(X) <- dt*grad
  } # for(ts in proj.steps)

  return(X.out)
}
