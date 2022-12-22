rm(list = ls())
library(dplyr)
library(ggplot2)
library(gridExtra)

source("coralsizedistfuncs.R")

#load data
oneyeardf <- read.csv("oneyeardf.csv",row.names=1)
oneyeardf <- oneyeardf[is.na(oneyeardf$ROI.LabelCode), ] # to remove corals that were out of frame 
axisscores <- read.csv("axisscores.csv", row.names=1)
axisscores <- axisscores[order(axisscores$PC1), ]

sites <- unique(oneyeardf$Site)
nsites <- length(sites)

PLB.bMLE.site.b <- numeric(nsites)
MSPLB.bMLE.site.b <- numeric(nsites)

siteb_plot <- list()

w <- 3648/35
v <- 2736/35

for(i in 1:nsites){
  s <- sites[i]
  sitedata <- oneyeardf %>% filter(Site == s)
  siteinput <- set.params(sitedata$Area)
  bML <- mle_b(Site == s, x = siteinput$Area, log_x = sitedata$log.Area, sum_log_x = siteinput$sum.log.Area,
                         x_min = siteinput$min.Area, x_max = siteinput$max.Area)
  PLB.bMLE.site.b[i] <- bML[[1]] 
  #PLB.minLL.site.b <- PLB.return.site[[2]]
  #PLB.minNegLL.site.b <- PLB.minLL.site.b$minimum
  MSPLB.bMLE.site.b[i] <- estimatebMSBPL(x = sitedata$Area, w = w, v = v)$minimum
  
  #x <- sitedata$Area
  sitex.PLB = seq(min(sitedata$Area), max(sitedata$Area), length = 1000)
  sitey.PLB = (1 - pPLB(x = sitex.PLB, b = PLB.bMLE.site.b[i], xmin = min(sitex.PLB),
                       xmax = max(sitex.PLB))) * length(sitedata$Area)
  siteb_plot[[i]] <- ggplot() +
    geom_point(aes(x = (sort(sitedata$Area, decreasing=TRUE)), y = (1:length(sitedata$Area))),
               color = "cadetblue", size = 2, alpha = 0.3) +
    xlab(expression(paste("Colony area, ", italic("x"), ~(cm^2)))) +
    ylab(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    "))) +
    scale_y_continuous(trans = 'log10', breaks = c(1,10,100), # currently capping at 100 because it is plotting Woolgoolga only; annoyingly the y axis doesnt scale but this will not be a problem when we plot the actual data and fit for each site - we will need to re-add 500 and 3000
                       limits = c(0.25, max(table(oneyeardf$Site)))) +
    scale_x_continuous(trans = 'log10', breaks = c(0,1,5,10,100,1000,10000),
                       limits = range(oneyeardf$Area))+
    geom_line(aes(x = sitex.PLB, y = sitey.PLB), col = 'black', lwd = 1) +
    annotate("text", x = 5, y = 10, label = s) +
    annotate("text", x = 5, y = 3, label = bquote(paste(italic(b)[PLB]==.(round(PLB.bMLE.site.b[i],2))))) +
    annotate("text", x = 5, y = 1, label = bquote(paste(italic(b)[MSBPL]==.(round(MSPLB.bMLE.site.b[i],2))))) +
    theme_classic()
}
do.call(grid.arrange, siteb_plot)
