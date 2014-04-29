############################
####  Projection model  ####
############################

NG <- 2   # number of sexes
AG <- 17  # number of age groups
DS <- 8   # number of disease stages
TS <- 4   # ART treatment duration stages
OUT.STEPS <- length(1970:2012)

dyn.load("rlib.so")
fnSpectrum <- function(iota, rVec){
  return(.C("fnSpectrumR", iota, rVec, OUT.STEPS, array(0, c(OUT.STEPS,  NG, AG, DS, TS)))[[4]])
}
