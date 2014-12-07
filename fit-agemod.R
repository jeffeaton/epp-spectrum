##################################################################
####                                                          ####
####  This script fits EPP and age-structured Spectrum model  ####
####  to test different assumptions about incorporating age-  ####
####  specific fertility and HIV sub-fertility into the data  ####
####  model for ANC data. R-spline model parameters are via   ####
####  maximum likelihood for 9 Southern and Eastern African   ####
####  countries.                                              ####
####                                                          ####
####  Results generated for report to UNAIDS Reference Group  ####
####  on Estiamtes, Modelling, and Projections on 5 December  ####
####  2014.                                                   ####
####                                                          ####
##################################################################


source("~/Documents/Code/R/read-epp-spectrum/read-spectrum-files.R")  # https://github.com/jeffeaton/read-epp-spectrum
source("~/Documents/Code/R/read-epp-spectrum/read-epp-files.R")

setwd("~/Documents/Code/R/anclik/")
source("anclik.R")

setwd("~/Documents/Research/age-specific-incidence/Spectrum/epp-spectrum/")

source("R/epp.R")
source("R/spectrum.R")
source("R/generics.R")

################################
####                        ####
####  HH survey likelihood  ####
####                        ####
################################

library(mvtnorm)

fnCreateSpectrumSubpops <- function(specfp, epp.input, epp.subpops, epp.data){

  ## This function distributes national number on ART to subnational
  ## regions.
  ## It takes a national specfp and epp.subpops and returns a specfp
  ## for each region.
  ## At present it simple (1) scales the init.pop for each region
  ## based on the 2010 population distribution into EPP, and
  ## (2) attributes ART according to average of regional prevalences
  ## as in EPP.
  ## It assumes the same fertility and mortality for each region
  ## as the national rates and no urbanisation / migration.

  ## Raise a warning if sum of subpops is more than 1% different from total population
  if(any(abs(1-rowSums(sapply(epp.subpops$subpops, "[[", "pop15to49")) / epp.input$epp.pop$pop15to49) > 0.01))
    warning("Sum of subpopulations does not equal total population")

  ancprev.means <- sapply(lapply(epp.data, "[[", "anc.prev"), mean, na.rm=TRUE)
  subpop.dist <- prop.table(sapply(epp.subpops$subpops, "[[", "pop15to49")[epp.subpops$total$year == 2010,])  # population distribution in 2010 (Dan's code)
  art.dist <- prop.table(subpop.dist * ancprev.means)

  specfp.subpop <- list()

  for(subpop in names(epp.subpops$subpops)){

    specfp.subpop[[subpop]] <- specfp
    specfp.subpop[[subpop]]$init.pop <- specfp$init.pop*subpop.dist[subpop]
    specfp.subpop[[subpop]]$artnum.15plus.ts <- specfp$artnum.15plus.ts*art.dist[subpop]

  }

  return(specfp.subpop)
}

fnPrepareHHSLikData <- function(hhs, anchor.year = 1970L){
  hhs$W.hhs <- qnorm(hhs$prev)
  hhs$v.hhs <- 2*pi*exp(hhs$W.hhs^2)*hhs$se^2  ## !! CHECK THIS (delta method)
  hhs$sd.W.hhs <- sqrt(hhs$v.hhs)
  hhs$idx <- hhs$year - (anchor.year - 1)

  hhslik.dat <- subset(hhs, used, c("W.hhs", "sd.W.hhs", "idx"))
  return(hhslik.dat)
}

fnCreateLikDat <- function(epp.data){

  return(list(anclik.dat = fnPrepareANCLikelihoodData(epp.data$anc.prev, epp.data$anc.n),
              hhslik.dat = fnPrepareHHSLikData(epp.data$hhs)))
}

fnHHSll <- function(qM, hhslik.dat){
  return(sum(dnorm(hhslik.dat$W.hhs, qM[hhslik.dat$idx], hhslik.dat$sd.W.hhs, log=TRUE)))
}

fnPregPrev.spec <- function(mod, fp){
  ## TODO: move this to spectrum.R file
  births.by.age <- t(fp$asfr.ts[,fp$proj.steps %% 1 == 0.5]) * rowSums(mod[,f.idx, fert.idx,,],,2)
  hivn.by.age <- rowSums(mod[,f.idx, fert.idx, 1,],,2)
  weighted.hivp.by.age <- sweep(rowSums(sweep(rowSums(mod[,f.idx, fert.idx, -1,1:3],,3), 3, fp$stage.fertrat, "*"),,2) + fp$art.fertrat * rowSums(mod[,f.idx, fert.idx,-1,4],,2), 2, fp$age.fertrat, "*")
  frac.hivp.moth <- weighted.hivp.by.age/(weighted.hivp.by.age+hivn.by.age)
  return(rowSums(births.by.age * frac.hivp.moth) / rowSums(births.by.age))
}


incid.spec <- function(mod, param, fp, age.idx=age15to49.idx, sex.idx=c(m.idx, f.idx)){
  param$rvec[fp$proj.steps %% 1 == 0.5] * (rowSums(mod[, , age15to49.idx, -1, 1]) + fp$relinfectART * rowSums(mod[, , age15to49.idx, -1, -1])) / rowSums(mod[,, age15to49.idx,1,])
}


ll <- function(theta, fp, likdat, pregprev=FALSE){
  theta.last <<- theta
  param <- fnCreateParam(theta, fp)

  if(inherits(fp, "specfp"))  ## TODO: revise to use generic functions
    mod <- fnSpectrum(param, fp)
  else
    mod <- fnEPP(param, fp)

  qM.all <- qnorm(prev(mod))
  qM.preg <- if(pregprev) qnorm(fnPregPrev(mod, fp)) else qM.all

  if(any(is.na(qM.all)) || any(qM.all[-1] == -Inf))
    return(-Inf)

  ll.anc <- log(fnANClik(qM.preg+param$ancbias, likdat$anclik.dat))
  ll.hhs <- fnHHSll(qM.all, likdat$hhslik.dat)
  return(ll.anc+ll.hhs)
}

fit.mod <- function(theta.init, fp, likdat, pregprev=FALSE, fp.update = list(), method="BFGS", maxit = 1000){
  fp <- do.call(update.specfp, c(list(fp), fp.update))
  fit <- optim(theta.init, ll, fp=fp, likdat=likdat, pregprev=pregprev, method=method, control=list(fnscale=-1, trace=4, REPORT=10, maxit = maxit))
  fit$theta.init <- theta.init
  fit$fp <- fp
  fit$param <- fnCreateParam(fit$par, fp)
  if(inherits(fp, "specfp"))
     fit$mod <- fnSpectrum(fit$param, fp)
  else
    fit$mod <- fnEPP(fit$param, fp)
  fit$prev <- prev(fit$mod)
  fit$incid <- incid(fit$mod, fit$param, fp)
  return(fit)
}


prepare.fit <- function(spec.path, upd.path, proj.end=2013.5){

  ## spectrum
  demp <- read.demog.param(upd.path)
  projp <- read.hivproj.param.5.0(paste(spec.path, ".DP", sep=""))

  specfp <- CreateSpectrumFixpar(projp, demp, proj.end = proj.end)

  ## epp
  eppd <- read.epp.data(paste(spec.path, ".xml", sep=""))
  epp.subp <- read.epp.subpops(paste(spec.path, ".xml", sep=""))
  epp.input <- read.epp.input(spec.path)

  epp.subp.input <- fnCreateEPPSubpops(epp.input, epp.subp, eppd)

  ## output
  val <- setNames(vector("list", length(eppd)), names(eppd))

  set.list.attr <- function(obj, attrib, value.lst)
    mapply(function(set, value){ attributes(set)[[attrib]] <- value; set}, obj, value.lst)

  val <- set.list.attr(val, "eppd", eppd)
  val <- set.list.attr(val, "likdat", lapply(eppd, fnCreateLikDat))
  val <- set.list.attr(val, "specfp", fnCreateSpectrumSubpops(specfp, epp.input, epp.subp, eppd))
  val <- set.list.attr(val, "eppfp", lapply(epp.subp.input, fnCreateEPPFixPar, proj.end = proj.end))
  val <- set.list.attr(val, "country", attr(eppd, "country"))

  return(val)
}


####  dev  ####

bw.upd <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/unpop/Botswana_72.upd"
ke.upd <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/unpop/Kenya_404.upd"
ls.upd <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/unpop/Lesotho_426.upd"
mw.upd <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/unpop/Malawi_454.upd"
tz.upd <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/unpop/United Republic of Tanzania_834.upd"
ug.upd <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/unpop/Uganda_800.upd"
za.upd <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/unpop/South Africa_710.upd"
zm.upd <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/unpop/Zambia_894.upd"
zw.upd <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/unpop/Zimbabwe_716.upd"

bw.path <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/Botswana 2014/Botswana 2014_Nat 19_06_14-c"
ke.path <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/Kenya 2014/Kenya 2014 May 5_c"
ls.path <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/Lesotho 2014/Lesotho May 12 2014-c"
mw.path <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/Malawi 2014/Malawi_2014_25June2014-c"
tz.path <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/UR\ Tanzania\ 2014/Tanzania_Mainland_Reviewed PMTCTApril 2014-c"
ug.path <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/Uganda 2014/Uganda RSPLINE 17th May 2014-c"
za.path <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/South Africa 2014/SouthAfrica_p_2014June10-H-c"
zm.path <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/Zambia 2014/Zambia 2014_Nat 31march 2014-mod-c"
zw.path <- "~/Documents/Data/Spectrum files/2014, final (downloaded 8 October 2014)/Zimbabwe 2014/Zimbabwe_23_April_2014_Finalupd-c"

bw.out <- prepare.fit(bw.path, bw.upd)
ke.out <- prepare.fit(ke.path, ke.upd, proj.end=2012.5)
ls.out <- prepare.fit(ls.path, ls.upd)
mw.out <- prepare.fit(mw.path, mw.upd, proj.end=2010.5)
tz.out <- prepare.fit(tz.path, tz.upd)
ug.out <- prepare.fit(ug.path, ug.upd, proj.end=2012.5)
za.out <- prepare.fit(za.path, za.upd)
zm.out <- prepare.fit(zm.path, zm.upd, proj.end=2011.5)
zw.out <- prepare.fit(zw.path, zw.upd)




####  Define fixed parameters for fitting  ####
base.fp <- list()
preg.fp <- list("cd4stage.weight" = c(1.3, 0.6, 0.1, 0.1, 0.0, 0.0, 0.0),  ## for EPP stage weighted prev
                "art1yr.weight" = 0.3)
asfr.fp <- list("age.fertrat" = rep(1, 7), "stage.fertrat" = rep(1, 7), "art.fertrat" = 1)
age.fp <- list("age.fertrat" = c(1.20, 0.76, 0.71, 0.65, 0.59, 0.53, 0.47), "stage.fertrat" = rep(1, 7), "art.fertrat" = 1)
stage.fp <- list("age.fertrat" = c(1.59, 1.06, 1.06, 1.00, 0.98, 0.91, 0.86),
                 "stage.fertrat" = c(1.25, 0.6, 0.4, 0.4, 0.3, 0.3, 0.3),
                 "art.fertrat" = 0.6)


## Botswana
theta.init <- c(-13.3952695, -7.3103391, -0.2410202, -0.5928517, -2.3756853, -3.3433023, 3.8777991, -3.1549190, 0.2638390)
bw.epp.fit <- lapply(bw.out, function(obj) fit.mod(theta.init, attr(obj, "eppfp"), attr(obj, "likdat"), fp.update=base.fp))
bw.preg.fit <- lapply(bw.out, function(obj) fit.mod(theta.init, attr(obj, "eppfp"), attr(obj, "likdat"), pregprev=TRUE, preg.fp))
bw.spec.fit <- lapply(bw.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), fp.update=base.fp))
bw.asfr.fit <- lapply(bw.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, asfr.fp))
bw.age.fit <- lapply(bw.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, age.fp))
bw.stage.fit <- lapply(bw.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, stage.fp))

## Kenya
theta.init <- c(-2.1796121, -1.4746154, -1.1848374, -0.1372349, -4.1435845, -1.4896543, -0.7685019, -4.6096204, 0)
ke.epp.fit <- lapply(ke.out, function(obj) fit.mod(theta.init, attr(obj, "eppfp"), attr(obj, "likdat"), fp.update=base.fp))

theta.init <- c(-4.6937959, 5.8595847, -12.4034238, 3.3067876, -5.9121974, 1.7473941, -10.2654680, -0.6503480, 0)
ke.preg.fit <- lapply(ke.out, function(obj) fit.mod(theta.init, attr(obj, "eppfp"), attr(obj, "likdat"), pregprev=TRUE, preg.fp, maxit=1000))

theta.init <- c(-2.1796121, -1.4746154, -1.1848374, -0.1372349, -4.1435845, -1.4896543, -0.7685019, -4.6096204, -0.1392529)
ke.spec.fit <- lapply(ke.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), fp.update=base.fp))
ke.asfr.fit <- lapply(ke.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, asfr.fp))
ke.age.fit <- lapply(ke.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, age.fp))
ke.stage.fit <- lapply(ke.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, stage.fp))

## Lesotho
theta.init <- c(rep(log(0.1), 7), log(0.01), 0)
ls.epp.fit <- lapply(ls.out, function(obj) fit.mod(theta.init, attr(obj, "eppfp"), attr(obj, "likdat"), fp.update=base.fp))
ls.preg.fit <- lapply(ls.out, function(obj) fit.mod(theta.init, attr(obj, "eppfp"), attr(obj, "likdat"), pregprev=TRUE, preg.fp))
ls.spec.fit <- lapply(ls.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), fp.update=base.fp))
ls.asfr.fit <- lapply(ls.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, asfr.fp))
ls.age.fit <- lapply(ls.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, age.fp))
ls.stage.fit <- lapply(ls.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, stage.fp))


## Malawi
theta.init <- list(c(-11.676080, -1.9100643, -0.008728004, -0.3291545, -3.080063, -1.035394, -18.54289, -10.442861, 0.3961385),
                   c(-11.676080, -1.9100643, -0.008728004, -0.3291545, -3.080063, -1.035394, -18.54289, -10.442861, 0.3961385),
                   c(-11.676080, -1.9100643, -0.008728004, -0.3291545, -3.080063, -1.035394, -18.54289, -10.442861, 0))
mw.epp.fit <- mapply(function(obj, theta.init) fit.mod(theta.init, attr(obj, "eppfp"), attr(obj, "likdat"), fp.update=base.fp), mw.out, theta.init, SIMPLIFY=FALSE)

theta.init <- list(c(-11.676080, -1.9100643, -0.008728004, -0.3291545, -3.080063, -1.035394, -18.54289, -10.442861, 0.3961385),
                   c(-2.1480645, -1.0021341, -1.0299901, -0.9127827, -2.9704270, -2.1800461, -1.4716459, -5.7241991, 0.5099688),
                   c(-11.676080, -1.9100643, -0.008728004, -0.3291545, -3.080063, -1.035394, -18.54289, -10.442861, 0))
mw.preg.fit <- mapply(function(obj, theta.init) fit.mod(theta.init, attr(obj, "eppfp"), attr(obj, "likdat"), pregprev=TRUE, fp.update=preg.fp), mw.out, theta.init, SIMPLIFY=FALSE)

theta.init <- list(c(-11.676080, -1.9100643, -0.008728004, -0.3291545, -3.080063, -1.035394, -18.54289, -10.442861, 0.3961385),
                   c(-11.676080, -1.9100643, -0.008728004, -0.3291545, -3.080063, -1.035394, -18.54289, -10.442861, 0.3961385),
                   c(-11.676080, -1.9100643, -0.008728004, -0.3291545, -3.080063, -1.035394, -18.54289, -10.442861, 0))
mw.spec.fit <- mapply(function(obj, theta.init) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), fp.update=base.fp), mw.out, theta.init, SIMPLIFY=FALSE)
mw.asfr.fit <- mapply(function(obj, theta.init) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, fp.update=asfr.fp), mw.out, theta.init, SIMPLIFY=FALSE)
mw.age.fit <- mapply(function(obj, theta.init) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, fp.update=age.fp), mw.out, theta.init, SIMPLIFY=FALSE)
mw.stage.fit <- mapply(function(obj, theta.init) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, fp.update=stage.fp), mw.out, theta.init, SIMPLIFY=FALSE)


## Tanzania
theta.init <- c(rep(log(0.15), 7), log(0.01), 0)
tz.epp.fit <- lapply(tz.out, function(obj) fit.mod(theta.init, attr(obj, "eppfp"), attr(obj, "likdat"), fp.update=base.fp))
tz.preg.fit <- lapply(tz.out, function(obj) fit.mod(theta.init, attr(obj, "eppfp"), attr(obj, "likdat"), pregprev=TRUE, preg.fp))
tz.spec.fit <- lapply(tz.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), fp.update=base.fp))
tz.asfr.fit <- lapply(tz.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, asfr.fp))
tz.age.fit <- lapply(tz.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, age.fp))
tz.stage.fit <- lapply(tz.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, stage.fp))

## Uganda
## theta.init <- c(0.1864551, -0.4043610, -0.5702416, -2.4875666, -2.5818859, -1.6024154, -3.6296576, -5.4023893, -0.1765377)

theta.init <- c(rep(log(0.1), 7), log(0.01), 0)
ug.epp.fit <- lapply(ug.out, function(obj) fit.mod(theta.init, attr(obj, "eppfp"), attr(obj, "likdat"), fp.update=base.fp))
ug.preg.fit <- lapply(ug.out, function(obj) fit.mod(theta.init, attr(obj, "eppfp"), attr(obj, "likdat"), pregprev=TRUE, preg.fp))
ug.spec.fit <- lapply(ug.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), fp.update=base.fp))
ug.asfr.fit <- lapply(ug.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, asfr.fp))
ug.age.fit <- lapply(ug.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, age.fp))
ug.stage.fit <- lapply(ug.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, stage.fp))


## South Africa
theta.init <- c(-9.3611031, -1.4070114, -1.2906695, 1.5985234, -3.3670246, 0.2470786, -15.6382398, -15.3049634, 0.3)
za.epp.fit <- lapply(za.out, function(obj) fit.mod(theta.init, attr(obj, "eppfp"), attr(obj, "likdat"), fp.update=base.fp))
za.preg.fit <- lapply(za.out, function(obj) fit.mod(theta.init, attr(obj, "eppfp"), attr(obj, "likdat"), pregprev=TRUE, preg.fp))
za.spec.fit <- lapply(za.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), fp.update=base.fp))
za.asfr.fit <- lapply(za.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, asfr.fp))
za.age.fit <- lapply(za.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, age.fp))
za.stage.fit <- lapply(za.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, stage.fp))

## Zambia
theta.init <- c(rep(log(0.15), 7), log(0.01), 0.4)
zm.epp.fit <- lapply(zm.out, function(obj) fit.mod(theta.init, attr(obj, "eppfp"), attr(obj, "likdat"), fp.update=base.fp))
zm.preg.fit <- lapply(zm.out, function(obj) fit.mod(theta.init, attr(obj, "eppfp"), attr(obj, "likdat"), pregprev=TRUE, preg.fp, maxit=10000))
zm.spec.fit <- lapply(zm.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), fp.update=base.fp))
zm.asfr.fit <- lapply(zm.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, asfr.fp))
zm.age.fit <- lapply(zm.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, age.fp))
zm.stage.fit <- lapply(zm.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, stage.fp))

## Zimbabwe
theta.init <- c(2.97218049, 8.37853230, -9.57620938, 1.84941942, -4.89303137, 1.51283083, -21.46272766, -43.39017534, 0.03120593)
zw.epp.fit <- lapply(zw.out, function(obj) fit.mod(theta.init, attr(obj, "eppfp"), attr(obj, "likdat"), fp.update=base.fp))
zw.preg.fit <- lapply(zw.out, function(obj) fit.mod(theta.init, attr(obj, "eppfp"), attr(obj, "likdat"), pregprev=TRUE, preg.fp))
zw.spec.fit <- lapply(zw.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), fp.update=base.fp))
zw.asfr.fit <- lapply(zw.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, asfr.fp))
zw.age.fit <- lapply(zw.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, age.fp))
zw.stage.fit <- lapply(zw.out, function(obj) fit.mod(theta.init, attr(obj, "specfp"), attr(obj, "likdat"), pregprev=TRUE, stage.fp))


save.image("workspace.RData")

compile.output <- function(out, country){

  attrs <- lapply(out, attributes)
  
  val <- mapply(list,
                "epp"   = get(paste(country, ".epp.fit", sep="")),
                "preg"  = get(paste(country, ".preg.fit", sep="")),
                "spec"  = get(paste(country, ".spec.fit", sep="")),
                "asfr"  = get(paste(country, ".asfr.fit", sep="")),
                "age"   = get(paste(country, ".age.fit", sep="")),
                "stage" = get(paste(country, ".stage.fit", sep="")),
                SIMPLIFY=FALSE)

  for(item in seq_along(val))
    attributes(val[[item]]) <- c(attributes(val[[item]]), attrs[[item]])

  return(val)
}

bw.out <- compile.output(bw.out, "bw")
ke.out <- compile.output(ke.out, "ke")
ls.out <- compile.output(ls.out, "ls")
mw.out <- compile.output(mw.out, "mw")
tz.out <- compile.output(tz.out, "tz")
ug.out <- compile.output(ug.out, "ug")
za.out <- compile.output(za.out, "za")
zm.out <- compile.output(zm.out, "zm")
zw.out <- compile.output(zw.out, "zw")

all.out <- c(bw.out, ke.out, ls.out, mw.out, tz.out, za.out, zm.out, zw.out)

for(region in seq_along(all.out))
  attr(all.out[[region]], "region") <- names(all.out)[region]

which(sapply(all.out, "attr", "country") == "United Republic of Tanzania")
attr(all.out[[10]], "country") <- "Tanzania"
attr(all.out[[11]], "country") <- "Tanzania"



###########################
####  Plot prevalence  ####
###########################

library(adegenet)

region <- all.out[[3]]
anchor.year <- 1970L

plot.prev <- function(region, anchor.year=1970L){

  likdat <- attr(region, "likdat")
  plot(NA, NA, type="l", ylim=c(0, 1.1*pnorm(max(unlist(likdat$anclik.dat$W.lst)) - region$epp$param$ancbias)),
       xlim=c(1990, 2013), lwd=2, xlab="Year", ylab="HIV prevalence", main=paste(attr(region, "country"), "-", attr(region, "region")))
  for(site in names(likdat$anclik.dat$W.lst)){
    points(likdat$anc$anc.idx.lst[[site]]+anchor.year-1, pnorm(likdat$anc$W.lst[[site]] - region$epp$param$ancbias), pch=17, col=transp("forestgreen", 0.3))
    lines(likdat$anc$anc.idx.lst[[site]]+anchor.year-1, pnorm(likdat$anc$W.lst[[site]] - region$epp$param$ancbias), col=transp("forestgreen", 0.3))
  }
  points(likdat$hhs$idx+anchor.year-1, pnorm(likdat$hhs$W.hhs), col=transp("blue", 0.85), lwd=2)
  segments(x0=likdat$hhs$idx+anchor.year-1,
           y0=pnorm(likdat$hhs$W.hhs - qnorm(0.975)*likdat$hhs$sd.W.hhs),
           y1=pnorm(likdat$hhs$W.hhs + qnorm(0.975)*likdat$hhs$sd.W.hhs),
           lwd=2, col=transp("blue", 0.85))
  xx <- 1969+1:length(region$epp$prev)
  lines(xx, region$epp$prev, col=1, lwd=1.5)
  lines(xx, region$preg$prev, col=2, lwd=1.5)
  lines(xx, region$spec$prev, col=1, lwd=1.5, lty=2)
  lines(xx, region$stage$prev, col=2, lwd=1.5, lty=2)
}


pdf("prevfit_2014-12-05.pdf", h=9, w=6, pointsize=7)

quartz(h=9, w=6, pointsize=7)
par(mfrow=c(4, 2), las=1, tcl=-0.25, mgp=c(1.8, 0.5, 0), mar=c(2, 3, 2, 0.5), cex=1.0)
lapply(all.out, plot.prev)
dev.off()


plot.incid <- function(region, anchor.year=1970L){

  max.incid <- max(region$epp$incid[-c(1:18)],
                   region$preg$incid[-c(1:18)],
                   region$spec$incid[-c(1:18)],
                   region$stage$incid[-c(1:18)])
        
  plot(NA, NA, type="l", ylim=c(0, 110*max.incid),
       xlim=c(1990, 2013), lwd=2, xlab="Year", ylab="HIV incidence", main=paste(attr(region, "country"), "-", attr(region, "region")))
  
  xx <- 1969+1:length(region$epp$prev)
  lines(xx, 100*region$epp$incid, col=1, lwd=1.5)
  lines(xx, 100*region$preg$incid, col=2, lwd=1.5)
  lines(xx, 100*region$spec$incid, col=1, lwd=1.5, lty=2)
  lines(xx, 100*region$stage$incid, col=2, lwd=1.5, lty=2)
    
}


pdf("incidfit_2014-12-05.pdf", h=9, w=6, pointsize=7)
par(mfrow=c(4, 2), las=1, tcl=-0.25, mgp=c(1.8, 0.5, 0), mar=c(2, 3, 2, 0.5), cex=1.0)
lapply(all.out, plot.incid)
dev.off()


plot.fit <- function(region, anchor.year=1970L){

  likdat <- attr(region, "likdat")
  plot(NA, NA, type="l", ylim=c(0, 1.1*pnorm(max(unlist(likdat$anclik.dat$W.lst)) - region$epp$param$ancbias)),
       xlim=c(1990, 2013), lwd=2, xlab="Year", ylab="HIV prevalence", main=paste(attr(region, "country"), "-", attr(region, "region")))
  for(site in names(likdat$anclik.dat$W.lst)){
    points(likdat$anc$anc.idx.lst[[site]]+anchor.year-1, pnorm(likdat$anc$W.lst[[site]] - region$epp$param$ancbias), pch=17, col=transp("forestgreen", 0.3))
    lines(likdat$anc$anc.idx.lst[[site]]+anchor.year-1, pnorm(likdat$anc$W.lst[[site]] - region$epp$param$ancbias), col=transp("forestgreen", 0.3))
  }
  points(likdat$hhs$idx+anchor.year-1, pnorm(likdat$hhs$W.hhs), col=transp("blue", 0.85), lwd=2)
  segments(x0=likdat$hhs$idx+anchor.year-1,
           y0=pnorm(likdat$hhs$W.hhs - qnorm(0.975)*likdat$hhs$sd.W.hhs),
           y1=pnorm(likdat$hhs$W.hhs + qnorm(0.975)*likdat$hhs$sd.W.hhs),
           lwd=2, col=transp("blue", 0.85))
  xx <- 1969+1:length(region$epp$prev)
  lines(xx, region$epp$prev, col=1, lwd=1.5)
  lines(xx, region$preg$prev, col=2, lwd=1.5)
  lines(xx, region$spec$prev, col=1, lwd=1.5, lty=2)
  lines(xx, region$stage$prev, col=2, lwd=1.5, lty=2)


  max.incid <- max(region$epp$incid[-c(1:18)],
                   region$preg$incid[-c(1:18)],
                   region$spec$incid[-c(1:18)],
                   region$stage$incid[-c(1:18)])
        
  plot(NA, NA, type="l", ylim=c(0, 110*max.incid),
       xlim=c(1990, 2013), lwd=2, xlab="Year", ylab="HIV incidence", main=paste(attr(region, "country"), "-", attr(region, "region")))
  
  xx <- 1969+1:length(region$epp$prev)
  lines(xx, 100*region$epp$incid, col=1, lwd=1.5)
  lines(xx, 100*region$preg$incid, col=2, lwd=1.5)
  lines(xx, 100*region$spec$incid, col=1, lwd=1.5, lty=2)
  lines(xx, 100*region$stage$incid, col=2, lwd=1.5, lty=2)
  
}

pdf("plotfit_2014-12-05.pdf", h=9, w=6, pointsize=7)
par(mfrow=c(4, 2), las=1, tcl=-0.25, mgp=c(1.8, 0.5, 0), mar=c(2, 3, 2, 0.5), cex=1.0)
lapply(all.out, plot.fit)
dev.off()





###################################
####  Comparison of model fit  ####
###################################

t(sapply(all.out, sapply, "[[", "value"))

table(apply(sapply(all.out, sapply, "[[", "value"), 2, which.max))  # which is best fit?
table(apply(sapply(all.out, sapply, "[[", "value")[1:2,], 2, which.max))  # which is better between epp and preg?
table(apply(sapply(all.out, sapply, "[[", "value")[c(1,3),], 2, which.max))  # which is better between epp and spec?
table(apply(sapply(all.out, sapply, "[[", "value")[c(1:2,6),], 2, which.max))  # is stage better than both epp and preg?


#########################
####  Model runtime  ####
#########################

library(microbenchmark)



all.out.bak <- all.out

timings <- NULL
for(reg in seq_along(all.out)){

  preg <- all.out[[reg]]$preg
  stage <- all.out[[reg]]$stage
  
  qM.preg <- qnorm(fnPregPrev(fnEPP(preg$param, preg$fp), preg$fp)) + preg$param$ancbias
  qM.stage <- qnorm(fnPregPrev(fnSpectrum(stage$param, stage$fp), stage$fp)) + stage$param$ancbias

  anclik.dat <- attributes(all.out[[reg]])$likdat$anclik.dat
  
  mb <- microbenchmark(fnEPP(preg$param, preg$fp),
                         fnSpectrum(stage$param, stage$fp),
                         fnANClik(qM.preg, anclik.dat),
                         times=1e3)

  timings <- rbind(timings, tapply(mb$time, mb$expr, sum)*1e-9)
  print(timings)
}

num.clinics <- sapply(lapply(lapply(lapply(all.out, attr, "likdat"), "[[", "anclik.dat"), "[[", "W.lst"), length)
num.prevobs <- sapply(lapply(lapply(lapply(lapply(all.out, attr, "likdat"), "[[", "anclik.dat"), "[[", "W.lst"), unlist), length)
anclik.time <- timings[,3]

num.clinics <- c(10, 14, 17, 23, 2, 15, 11, 18, 25, 96, 93, 4, 8, 12, 7, 9, 6, 7, 5, 6, 12, 12, 35, 41)
anclik.time <- c(6.2144487, 5.6004745, 5.7982736, 1.4907405, 0.4942408, 3.9054802, 2.1441295, 3.1546658, 3.9167895, 6.2066913, 6.2419659, 2.3461965, 5.4722075, 9.7643534,
                 4.7767987, 6.7112732, 4.2531070, 4.7203950, 1.0693049, 3.4726152, 4.3240463, 2.6973985, 2.7774486, 2.7322841)


