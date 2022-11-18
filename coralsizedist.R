rm(list = ls())
library(dplyr)
library(ggplot2)
library(svglite)

source("functionsfromcarvalho.R")

#load data
oneyeardf <- read.csv("oneyeardf.csv",row.names=1)
oneyeardf <- oneyeardf[is.na(oneyeardf$ROI.LabelCode), ] # to remove corals that were out of frame 
axisscores <- read.csv("axisscores.csv", row.names=1)
axisscores <- axisscores[order(axisscores$PC1), ]

sites <- unique(oneyeardf$Site)

for(s in sites){
  sitedata <- oneyeardf %>% filter(Site == s)
  siteinput <- set.params(sitedata$Area)
  PLB.return.site <- mle_b(Site == s, x = siteinput$Area, log_x = sitedata$log.Area, sum_log_x = siteinput$sum.log.Area,
                         x_min = siteinput$min.Area, x_max = siteinput$max.Area)
  PLB.bMLE.site.b <- PLB.return.site[[1]] 
  PLB.minLL.site.b <- PLB.return.site[[2]]
  PLB.minNegLL.site.b <- PLB.minLL.site.b$minimum
  x <- siteinput$Area
  sitex.PLB = seq(min(siteinput$Area), max(siteinput$Area), length=1000)
  sitey.PLB = (1 - pPLB(x = sitex.PLB, b = PLB.bMLE.site.b, xmin = min(sitex.PLB),
                       xmax = max(sitex.PLB))) * length(siteinput$Area)
  #spectra.text <- as.character(round(PLB.bMLE.site.b, 2))
  siteb_plot <- ggplot() +
    geom_point(aes(x = (sort(siteinput$Area, decreasing=TRUE)), y = (1:length(siteinput$Area))),
               color = "cadetblue", size = 2, alpha = 0.3) +
    xlab(expression(paste("Colony area, ", italic("x"), " (cm^2)"))) +
    ylab(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    "))) +
    scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
                       limits = c(0.25, max(table(oneyeardf$Site)))) +
    scale_x_continuous(trans = 'log10', breaks = c(0,1,5,10,100,1000,10000),
                       limits = range(oneyeardf$Area))+
    geom_line(aes(x = sitex.PLB, y = sitey.PLB), col = 'black', lwd = 1) +
    annotate("text", x = 5, y = 10, label = s) +
    annotate("text", x = 5, y = 3, label = expression(paste(italic("b = "), bquote.(PLB.bMLE.site.b))) +
    theme_classic()
}

# fix the b value
# work out how to store all the numbers
