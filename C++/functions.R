############################
####  Projection model  ####
############################

NG <- 2   # number of sexes
AG <- 17  # number of age groups
DS <- 8   # number of disease stages
TS <- 4   # ART treatment duration stages
OUT.STEPS <- length(1970:2012)

age15plus.idx <- 4:AG

system("rm C++/*.o")
system("R CMD SHLIB -lgsl -lgslcblas -lgfortran C++/rlib.cpp C++/model.cpp C++/states.cpp C++/incidence.cpp C++/mvndstpack.f")

dyn.load("C++/rlib.so")

fnSpectrum <- function(param){
  return(.Call("fnSpectrumR", param, OUT.STEPS))
}


create.rmat <- function(m.mean, m.sd, f.mean, f.sd, corr, margdist = "gamma", copula = "gaussian"){

  if(margdist != "gamma" || copula != "gaussian")
    stop("not implemented")

  return(.Call("createRmatR", m.mean, m.sd, f.mean, f.sd, corr))
}
