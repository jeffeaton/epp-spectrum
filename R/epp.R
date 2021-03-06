library(splines)

fnCreateEPPSubpops <- function(epp.input, epp.subpops, epp.data){

  ## Raise a warning if sum of subpops is more than 1% different from total population
  if(any(abs(1-rowSums(sapply(epp.subpops$subpops, "[[", "pop15to49")) / epp.input$epp.pop$pop15to49) > 0.01))
    warning("Sum of subpopulations does not equal total population")

  ## EPP implmentation from Dan Hogan's code distributes ART according to relative average ANC prevalence in each subpopulation
  ## A more representative approach might be to apportion this based on national survey prevalence in each subpopulation
  ancprev.means <- sapply(lapply(epp.data, "[[", "anc.prev"), mean, na.rm=TRUE)
  subpop.dist <- prop.table(sapply(epp.subpops$subpops, "[[", "pop15to49")[epp.subpops$total$year == 2010,])  # population distribution in 2010 (Dan's code)
  art.dist <- prop.table(subpop.dist * ancprev.means)

  epp.subpop.input <- list()

  for(subpop in names(epp.subpops$subpops)){

    epp.subpop.input[[subpop]] <- epp.input
    epp.subpop.input[[subpop]]$epp.pop <- epp.subpops$subpops[[subpop]]

    epp.art <- epp.input$epp.art
    epp.art$m.val[epp.art$m.isperc == "N"] <- epp.art$m.val[epp.art$m.isperc == "N"] * art.dist[subpop]
    epp.art$f.val[epp.art$f.isperc == "N"] <- epp.art$f.val[epp.art$f.isperc == "N"] * art.dist[subpop]

    epp.subpop.input[[subpop]]$epp.art <- epp.art
  }

  return(epp.subpop.input)
}


fnCreateEPPFixPar <- function(epp.input, proj.start = epp.input$start.year, proj.end = epp.input$stop.year, dt = 0.1,
                              tsEpidemicStart = proj.start, cd4stage.weights=c(1.3, 0.6, 0.1, 0.1, 0.0, 0.0, 0.0), art1yr.weight = 0.1){

  #########################
  ##  Population inputs  ##
  #########################

  epp.pop <- epp.input$epp.pop
  proj.steps <- seq(proj.start, proj.end, dt)
  epp.pop.ts <- data.frame(pop15to49 = approx(epp.pop$year, epp.pop$pop15to49, proj.steps)$y,
                           age15enter = approx(epp.pop$year, dt*epp.pop$pop15, proj.steps)$y,
                           age50exit = approx(epp.pop$year, dt*epp.pop$pop50, proj.steps)$y,
                           netmigr = approx(epp.pop$year, dt*epp.pop$netmigr, proj.steps)$y)

  proj.years <- floor(proj.steps)

  epp.pop.ts$netmigr <- epp.pop.ts$netmigr * with(subset(epp.pop, epp.pop$year %in% proj.years), rep(ifelse(netmigr != 0, netmigr * table(proj.years) * dt / tapply(epp.pop.ts$netmigr, floor(proj.steps), sum), 0), times=table(proj.years)))
  epp.pop.ts$age50rate <- (epp.pop.ts$age50exit/epp.pop.ts$pop15to49)/dt
  epp.pop.ts$mx <- c(1.0 - (epp.pop.ts$pop15to49[-1] - (epp.pop.ts$age15enter - epp.pop.ts$age50exit + epp.pop.ts$netmigr)[-length(proj.steps)]) / epp.pop.ts$pop15to49[-length(proj.steps)], NA) / dt # SLIGHTLY DIFFERENT FROM DAN'S CODE: take age 50 entrants and exit from current DT rather than next one.
  epp.pop.ts[length(proj.steps), "mx"] <- epp.pop.ts[length(proj.steps)-1, "mx"]


  ##################
  ##  ART inputs  ##
  ##################

  epp.art <- epp.input$epp.art

  ## number of persons who should be on ART at end of timestep
  artnum.ts <- with(subset(epp.art, m.isperc=="N"), approx(year+1-dt, (m.val+f.val)*(1.0-perc50plus/100), proj.steps, rule=2))$y  # offset by 1 year because number on ART are Dec 31
  ## !! Currently only dealing with numbers on ART, assuming no increase in coverage after last number

  epp.art$artelig.idx <- match(epp.art$cd4thresh, c(1, 2, 500, 350, 250, 200, 100, 50))
  artelig.idx.ts <- approx(epp.art$year, epp.art$artelig.idx, proj.steps, "constant", rule=2)$y

  cd4prog <- 1/colMeans(epp.input$cd4stage.dur[c(2:3,6:7),])
  cd4init <- 0.01*colMeans(epp.input$cd4initperc[c(2:3,6:7),])
  cd4artmort <- cbind(colMeans(epp.input$cd4mort[c(2:3,6:7),]),
                      colMeans(epp.input$artmort.less6mos[c(2:3,6:7),]),
                      colMeans(epp.input$artmort.6to12mos[c(2:3,6:7),]),
                      colMeans(epp.input$artmort.after1yr[c(2:3,6:7),]))

  relinfectART <- 1.0 - epp.input$infectreduc


  ###########################
  ##  r-spline parameters  ##
  ###########################

  numKnots <- 7
  proj.dur <- diff(range(proj.steps))
  rvec.knots <- seq(min(proj.steps) - 3*proj.dur/(numKnots-3), max(proj.steps) + 3*proj.dur/(numKnots-3), proj.dur/(numKnots-3))
  rvec.spldes <- splineDesign(rvec.knots, proj.steps)

  val <- list(proj.steps      = proj.steps,
              tsEpidemicStart = tsEpidemicStart,
              dt              = dt,
              epp.pop.ts      = epp.pop.ts,
              artnum.ts       = artnum.ts,
              artelig.idx.ts  = artelig.idx.ts,
              cd4prog         = cd4prog,
              cd4init         = cd4init,
              cd4artmort      = cd4artmort,
              relinfectART    = relinfectART,
              numKnots        = numKnots,
              rvec.spldes     = rvec.spldes)

  class(val) <- "eppfp"
  return(val)
}



############################################################
############################################################

dyn.load("src/epp.so")

"incr<-" <- function(x, value) { x + value } # increment operator, from Hmisc

DS <- 8   # number of disease stages
TS <- 4   # ART treatment duration stages


fnEPP <- function(param, fp, VERSION = "C"){


  if(VERSION != "R"){
    mod <- .Call("fnEPP", fp$epp.pop.ts, fp$proj.steps, fp$dt,
                 param$rvec, param$iota, fp$relinfectART, as.numeric(fp$tsEpidemicStart),
                 fp$cd4init, fp$cd4prog, fp$cd4artmort,
                 fp$artnum.ts, as.integer(fp$artelig.idx.ts))
    class(mod) <- "epp"
    return(mod)
  }

##################################################################################

  proj.steps <- fp$proj.steps
  epp.pop.ts <- fp$epp.pop.ts
  dt <- fp$dt

  ## initialize output
  Xout <- array(NA, c(sum(fp$proj.steps %% 1 == dt*floor((1/dt)/2)), DS, TS))

  ## initialize population
  X <- array(0, c(DS, TS))
  X[1,1] <- epp.pop.ts$pop15to49[1]


  for(ts in 1:length(proj.steps)){

    if(proj.steps[ts] %% 1 == dt*floor((1/dt)/2)) # store the result mid-year (or ts before mid-year)
      Xout[ceiling(ts*dt),,] <- X

    grad <- array(0, c(DS, TS))

    ## ageing and natural mortality
    incr(grad) <- -X * (epp.pop.ts$age50rate[ts] + epp.pop.ts$mx[ts])

    ## new entrants
    incr(grad[1,1]) <- epp.pop.ts$age15enter[ts] / dt  # convert to annual rate

    ## net migrants
    incr(grad) <- X/sum(X) * epp.pop.ts$netmigr[ts] / dt

    ## new infections
    incrate <- param$rvec[ts] * (sum(X[-1,1]) + fp$relinfectART * sum(X[-1,-1])) / sum(X) + param$iota*(proj.steps[ts] == fp$tsEpidemicStart)
    incr(grad[1,1]) <- -X[1,1] * incrate
    incr(grad[-1,1]) <- X[1,1] * incrate * fp$cd4init

    ## disease progression and mortality
    incr(grad[2:(DS-1),1]) <- -fp$cd4prog * X[2:(DS-1),1]    # remove cd4 stage progression (untreated)
    incr(grad[3:DS,1]) <- fp$cd4prog * X[2:(DS-1),1]         # add cd4 stage progressiogn (untreated)
    incr(grad[-1, 2:3]) <- -2 * X[-1, 2:3]                # remove ART duration progression (HARD CODED 6 months duration)
    incr(grad[-1, 3:4]) <- 2 * X[-1, 2:3]                 # add ART duration progression (HARD CODED 6 months duration)
    incr(grad[-1,]) <- - fp$cd4artmort * X[-1,]              # HIV mortality

    ## ART initiation
    if(fp$artnum.ts[ts] > 0){
      artnum.curr <- sum(X[-1,-1])
      artnum.anninits <- (fp$artnum.ts[ts] - artnum.curr)/dt - sum(grad[-1,-1]) # desired change rate minus current exits

      artelig <- X[fp$artelig.idx.ts[ts]:DS,1]
      expect.mort.weight <- fp$cd4artmort[fp$artelig.idx.ts[ts]:DS - 1, 1] / sum(artelig * fp$cd4artmort[fp$artelig.idx.ts[ts]:DS - 1, 1])
      artinit.weight <- (expect.mort.weight + 1/sum(artelig))/2  # average eligibility and expected mortality
      artinit.ann.ts <- pmin(artnum.anninits * artinit.weight * artelig,
                             artelig/dt,                                   # check that don't initiate more than the number eligible
                             artelig/dt + grad[fp$artelig.idx.ts[ts]:DS,1], na.rm=TRUE)

      incr(grad[fp$artelig.idx.ts[ts]:DS, 1]) <- -artinit.ann.ts
      incr(grad[fp$artelig.idx.ts[ts]:DS, 2]) <- artinit.ann.ts
    }

    ## do projection (euler integration) ##
    incr(X) <- dt*grad

  } # for(ts in proj.steps)

  class(Xout) <- "epp"
  return(Xout)
}

prev.epp <- function(mod){
  return(rowSums(mod[,-1,])/rowSums(mod))
}

fnPregPrev.epp <- function(mod, fp){

  pregweight.hivn <- mod[,1,1]
  pregweight.hivp <- rowSums(sweep(mod[,-1,1:3], 2, fp$cd4stage.weight, "*"))
  pregweight.art1yr <- rowSums(mod[,-1,4]) * fp$art1yr.weight

  return((pregweight.hivp+pregweight.art1yr)/(pregweight.hivn+pregweight.hivp+pregweight.art1yr))
}

incid.epp <- function(mod, param, fp){
  param$rvec[fp$proj.steps %% 1 == 0.5] * (rowSums(mod[,-1,1]) + fp$relinfectART * rowSums(mod[,-1,-1])) / rowSums(mod)
}

## fnCreateParam <- function(theta, fp){
##   return(list(rvec = exp(as.vector(fp$rvec.spldes %*% theta[1:fp$numKnots])),
##               iota = exp(theta[fp$numKnots+1]),
##               ancbias = theta[fp$numKnots+2]))
## }

fnARTCov <- function(mod){
  rowSums(mod[,-1,-1]) / rowSums(mod[,-1,])
}

fnPregARTCov <- function(mod, cd4stage.weights=c(1.3, 0.6, 0.1, 0.1, 0.0, 0.0, 0.0), art1yr.weight = 0.3){

  pregweight.noart <- rowSums(sweep(mod[,-1,1], 2, cd4stage.weights, "*"))
  pregweight.art.less1yr <- rowSums(sweep(mod[,-1,2:3], 2, cd4stage.weights, "*"))
  pregweight.art1yr <- rowSums(mod[,-1,4]) * art1yr.weight

  return((pregweight.art.less1yr+pregweight.art1yr)/(pregweight.art.less1yr+pregweight.art1yr+pregweight.noart))
}
