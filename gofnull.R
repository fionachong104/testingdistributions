rm(list = ls())
source("coralsizedistfuncs.R")

gof <- read.table(file = "gofall.csv", header = TRUE, sep = ",")

par(mfrow = c(2, 2))
plotFisherapprox(X2 = gof$MSlognormX2, df = gof$MSlognormdf, label = "a: minus-sampled lognormal")
plotFisherapprox(X2 = gof$lognormX2, df = gof$lognormdf, label = "b: lognormal")
plotFisherapprox(X2 = gof$MSBPLX2, df = gof$MSBPLdf, label = "c: minus-sampled bounded power law")
plotFisherapprox(X2 = gof$BPLX2, df = gof$BPLdf, label = "d: bounded power law")
     