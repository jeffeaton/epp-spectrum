setwd("..")

source("R/incidence-model.R")
source("C++/functions.R")

round(create.Rmat(35, 20, 25, 15, 0.6), 4)

round(create.rmat(35, 20, 25, 15, 0.6), 4)

round(create.Rmat(35, 20, 25, 15, 0.6), 4) == round(create.rmat(35, 20, 25, 15, 0.6), 4)

system.time(replicate(100, create.Rmat(35, 20, 25, 15, 0.6)))
system.time(replicate(1000, create.rmat(35, 20, 25, 15, 0.6)))
