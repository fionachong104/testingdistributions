rm(list = ls())
source("coralsizedistfuncs.R")

gof <- read.table(file = "gofall.csv", header = TRUE, sep = ",")

par(mfrow = c(2, 1))
plotFisherapprox(X2 = gof$lognormX2, df = gof$lognormdf, label = "a: lognormal")
plotFisherapprox(X2 = gof$BPLX2, df = gof$BPLdf, label = "b: bounded power law")
     