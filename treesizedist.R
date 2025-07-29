rm(list = ls())
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)

source("coralsizedistfuncs.R")

#load data
tree <- read.csv("C:/Users/tvkx991/OneDrive - University of Leeds/PhD/testingdistributions/testingdistributions/Ausplots_large_tree_survey_data.csv")
tree <- tree %>% drop_na(aboveGroundLiveBiomass_kilograms)  # to remove dead trees
# axisscores <- read.csv("axisscores.csv", row.names=1)
# axisscores <- axisscores[order(axisscores$PC1), ]

sites <- unique(tree$siteName)
#sites <- sites[match(row.names(axisscores), sites)] #re-ordering based on increasing PC1 scores
nsites <- length(sites)

# create empty lists to store things in 
PLB.bMLE.site.b <- numeric(nsites)
siteb_plot <- list()
bplqq <- list()
AICdf <- data.frame(site = sites, llBPL = NA, lllognorm = NA, AICBPL = NA, AIClognorm = NA)


for(i in 1:nsites){
  s <- sites[i]
  sitedata <- tree %>% filter(siteName == s)
  siteinput <- set.fish.params(sitedata$aboveGroundLiveBiomass_kilograms)
  bML <- mle_b(Site == s, x = siteinput$biomass, sum_log_x = siteinput$sum.log.biomass,
               x_min = siteinput$min.biomass, x_max = siteinput$max.biomass)
  PLB.bMLE.site.b[i] <- bML[[1]] 
  thetalnorm <- estimatelognormal(x = sitedata$aboveGroundLiveBiomass_kilograms)
  x <- sitedata$aboveGroundLiveBiomass_kilograms
  sitex = seq(min(sitedata$aboveGroundLiveBiomass_kilograms), max(sitedata$aboveGroundLiveBiomass_kilograms), length = 10000)
  par(mfrow=c(1,2))
  bplqq[[i]] <- qqplot(FXinv(u = ppoints(siteinput$n), b = PLB.bMLE.site.b[i], xmin = min(sitex),
                             xmax = max(sitex)), x , xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(A)", sites[i], ": Power law Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
  #qqline(x, distribution = function(p){
  #  FXinv(p, b = PLB.bMLE.site.b[i], xmin = min(sitex), xmax = max(sitex))
  #}, untf=T)
  qqplot(qlnorm(p = ppoints(siteinput$n), meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog), x, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(B)", sites[i], ": Log-normal Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
  #  qqline(x, distribution = function(p){qlnorm(p, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog)}, untf=T)
  #rank plot   
  sitey.PLB <- (1 - pPLB(x = sitex, b = PLB.bMLE.site.b[i], xmin = min(sitex),
                         xmax = max(sitex))) * length(sitedata$aboveGroundLiveBiomass_kilograms)
  sitey.lnorm <- plnorm(q = sitex, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog, lower.tail = FALSE, log.p = FALSE) * length(sitedata$aboveGroundLiveBiomass_kilograms)
  siteb_plot[[i]] <- ggplot() +
    geom_point(aes_(x = (sort(sitedata$aboveGroundLiveBiomass_kilograms, decreasing=TRUE)), y = (1:length(sitedata$aboveGroundLiveBiomass_kilograms))),
               color = "#666666", size = 2, alpha = 0.3) +
    scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000), 
                       limits = c(0.25, max(table(tree$siteName)))) +
    scale_x_continuous(trans = 'log10',#breaks = c(-10,0,1,10,100),
                       limits = range(tree$aboveGroundLiveBiomass_kilograms))+
    geom_line(aes_(x = sitex, y = sitey.PLB), col = 'black', lwd = 1) +
    geom_line(aes_(x = sitex, y = sitey.lnorm), col = '#1B9E77', lwd = 1) +
    labs(tag = LETTERS[i]) +
    annotate("text", x = 0.002, y = 10, label = s) +
    annotate("text", x = 100, y = 100, label = bquote(paste(italic(b)[PLB]==.(round(PLB.bMLE.site.b[i],2))))) +
    annotate("text", x = 100, y =1000, label = bquote(n == .(length(sitedata$aboveGroundLiveBiomass_kilograms)))) +
    theme_classic() + 
    theme(axis.title = element_blank())
  AICdf$lllognorm[i] <- lnormAIC(x)$lllognorm
  AICdf$AIClognorm[i] <- lnormAIC(x)$AIClognorm
  AICdf$llBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.biomass, xmax = siteinput$max.biomass, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$llBPL
  AICdf$AICBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.biomass, xmax = siteinput$max.biomass, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$AICBPL
}

leftlabel <- grid::textGrob(expression(paste("Number of tree with sizes", " ">=" ", italic("x"), "    ")), rot = 90)
#bottomlabel <- grid::textGrob(expression(paste("tree biomass, ", italic("x"), ~(kg))))
bottomlabel <- grid::textGrob(expression(paste("tree biomass, ", italic("x"), ~(kg))))

# siteb_plot <- grid.arrange(grobs = siteb_plot, ncol = 2,
#                            left = leftlabel,
#                            bottom = bottomlabel)

marrangeGrob(grobs= siteb_plot, nrow=4, ncol=4)


write.csv(AICdf,'tree_AIC.csv')

# # # saves the temp images from the plotting envrionment
# plots.dir.path <- list.files(tempdir(), pattern = "rs-graphics", full.names = TRUE)
# plots.png.paths <- list.files(plots.dir.path, pattern = ".png", full.names = TRUE)
# file.copy(from = plots.png.paths, to = "C:/Users/624225/OneDrive - hull.ac.uk/_BoxData/PhD/testingdistributions/qqplots/fish_logged")
