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
  MSPLB.bMLE.site.b[i] <- estimatebMSBPL(x = sitedata$Area, w = w, v = v)$minimum
  #x <- sitedata$Area
  sitex = seq(min(sitedata$Area), max(sitedata$Area), length = 1000)
  sitey.PLB = (1 - pPLB(x = sitex, b = PLB.bMLE.site.b[i], xmin = min(sitex),
                       xmax = max(sitex))) * length(sitedata$Area)
  sitey.MSBPL = (1 - FMSBPL(x = sitex, b = MSPLB.bMLE.site.b[i], xmin = min(sitex),
                            w = w, v = v)) * length(sitedata$Area)
  siteb_plot[[i]] <- ggplot() +
    geom_point(aes_(x = (sort(sitedata$Area, decreasing=TRUE)), y = (1:length(sitedata$Area))),
               color = "cadetblue", size = 2, alpha = 0.3) +
    scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000), 
                       limits = c(0.25, max(table(oneyeardf$Site)))) +
    scale_x_continuous(trans = 'log10', breaks = c(0,1,5,10,100,1000,10000),
                       limits = range(oneyeardf$Area))+
    geom_line(aes_(x = sitex, y = sitey.PLB), col = 'black', lwd = 1) +
    geom_line(aes_(x = sitex, y = sitey.MSBPL), col = 'red', lwd = 1) +
    annotate("text", x = 10, y = 10, label = s) +
    annotate("text", x = 10, y = 3, label = bquote(paste(italic(b)[PLB]==.(round(PLB.bMLE.site.b[i],2))))) +
    annotate("text", x = 10, y = 1, label = bquote(paste(italic(b)[MSBPL]==.(round(MSPLB.bMLE.site.b[i],2))))) +
    theme_classic() + 
    theme(axis.title = element_blank())
}

leftlabel <- grid::textGrob(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    ")), rot = 90)
bottomlabel <- grid::textGrob(expression(paste("Colony area, ", italic("x"), ~(cm^2))))

grid.arrange(grobs = siteb_plot, ncol = 4, 
             left = leftlabel,
             bottom = bottomlabel)

