library(gdata)
fert.dist <- read.xls("../south-africa/south-africa-demproj.xlsx", "ASFR", row.names=1)
tfr <- as.numeric(read.xls("../south-africa/south-africa-demproj.xlsx", "TFR", row.names=1))
asfr <- mapply("*", fert.dist[1:7,]/500, as.list(tfr))
fert.rat <- read.xls("../south-africa/south-africa-nathist.xlsx", "fertility ratio", row.names=1)$Ratio
rm(fert.dist, tfr)

source("functions.R")


## Need to finish writing likelihood function ##
