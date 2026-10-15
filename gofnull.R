rm(list = ls())
source("coralsizedistfuncs.R")

gof <- read.table(file = "gofall.csv", header = TRUE, sep = ",")

coral <- gof[1:20,]

par(mfrow = c(2, 1))
plotFisherapprox(X2 = coral$lognormX2, df = coral$lognormdf, label = "a: lognormal", obsperbin = 5)
plotFisherapprox(X2 = coral$BPLX2, df = coral$BPLdf, label = "b: bounded power law", obsperbin = 5)
     
plotFisherapprox(X2 = gof$lognormX2, df = gof$lognormdf, label = "a: lognormal")
plotFisherapprox(X2 = gof$BPLX2, df = gof$BPLdf, label = "b: bounded power law")
