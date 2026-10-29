rm(list = ls())
source("coralsizedistfuncs.R")

gof <- read.table(file = "gofall.csv", header = TRUE, sep = ",")

coral <- gof[1:20,]
fish <- gof[21:43,]

par(mfrow = c(2, 1))
plotFisherapprox(X2 = coral$lognormX2, df = coral$lognormdf, label = "A: Log-normal", obsperbin = 5)
plotFisherapprox(X2 = coral$BPLX2, df = coral$BPLdf, label = "B: Bounded power law", obsperbin = 5)
     
plotFisherapprox(X2 = gof$lognormX2, df = gof$lognormdf, label = "A: Log-normal")
plotFisherapprox(X2 = gof$BPLX2, df = gof$BPLdf, label = "B: Bounded power law")

plotFisherapprox(X2 = fish$lognormX2, df = fish$lognormdf, label = "A: Log-normal", obsperbin = 5)
plotFisherapprox(X2 = fish$BPLX2, df = fish$BPLdf, label = "B: Bounded power law", obsperbin = 5)
