rm(list = ls())
source("coralsizedistfuncs.R")

gof <- read.table(file = "gofall.csv", header = TRUE, sep = ",")

plotFisherapprox(X2 = gof$MSlognormX2, df = gof$MSlognormdf)
     