rm(list = ls())
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)

source("coralsizedistfuncs.R")

#load data
bird <- read.csv("Weights.csv")
#tree <- tree %>% drop_na(aboveGroundLiveBiomass_kilograms)  # to remove dead trees
# axisscores <- read.csv("axisscores.csv", row.names=1)
# axisscores <- axisscores[order(axisscores$PC1), ]
bird$siteName = paste(bird$location, bird$Year, sep = "_")
sites <- unique(bird$siteName)
#sites <- sites[match(row.names(axisscores), sites)] #re-ordering based on increasing PC1 scores
nsites <- length(sites)

# create empty lists to store things in 
PLB.bMLE.site.b <- numeric(nsites)
siteb_plot <- list()
AICdf <- data.frame(Site = sites, Number = NA, llBPL = NA, lllognorm = NA, AICBPL = NA, AIClognorm = NA)
gof <- data.frame(site = sites, lognormX2 = NA, lognormdf = NA, lognormP = NA, BPLX2 = NA, BPLdf = NA, BPLP = NA)


for(i in 1:nsites){
  s <- sites[i]
  sitedata <- bird %>% filter(siteName == s)
  siteinput <- set.bird.params(sitedata$Weight)
  bML <- mle_b(Site == s, x = siteinput$biomass, sum_log_x = siteinput$sum.log.biomass,
               x_min = siteinput$min.biomass, x_max = siteinput$max.biomass)
  PLB.bMLE.site.b[i] <- bML[[1]]
  thetalnorm <- estimatelognormal(x = sitedata$Weight)

  #goodness-of-fit test for lognormal
  lngof <- lnormgof(x = sitedata$Weight, mu = thetalnorm$meanlog, sigma = thetalnorm$sdlog)
  gof$lognormX2[i] <- lngof$X2
  gof$lognormdf[i] <- lngof$df
  gof$lognormP[i] <- lngof$P

  #goodness-of-fit test for bounded power law
  bplgof <- BPLgof(x = sitedata$Weight, b = bML[[1]], xmin = min(sitedata$Weight), xmax = max(sitedata$Weight))
  gof$BPLX2[i] <- bplgof$X2
  gof$BPLdf[i] <- bplgof$df
  gof$BPLP[i] <- bplgof$P

  x <- sitedata$Weight
  sitex = seq(min(sitedata$Weight), max(sitedata$Weight), length = 10000)
    #rank plot   
  sitey.PLB <- (1 - pPLB(x = sitex, b = PLB.bMLE.site.b[i], xmin = min(sitex),
                         xmax = max(sitex))) * length(sitedata$Weight)
  sitey.lnorm <- plnorm(q = sitex, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog, lower.tail = FALSE, log.p = FALSE) * length(sitedata$Weight)
  siteb_plot[[i]] <- ggplot() +
    geom_point(aes_(x = (sort(sitedata$Weight, decreasing=TRUE)), y = (1:length(sitedata$Weight))),
               color = "#666666", size = 2, alpha = 0.3) +
    scale_y_continuous(trans = 'log10', breaks = c(1,10,100,350), 
                       limits = c(0.25, max(table(bird$siteName)))) +
    scale_x_continuous(trans = 'log10', breaks = c(10, 100, 500),
                       limits = range(bird$Weight))+
    geom_line(aes_(x = sitex, y = sitey.PLB), col = 'black', lwd = 1) +
    geom_line(aes_(x = sitex, y = sitey.lnorm), col = '#1B9E77', lwd = 1) +
    labs(tag = LETTERS[ ((i - 1) %% 4) + 1 ])+
     annotate("text", x = 10, y = 10, label = s) +
  annotate("text", x = 10, y = 5, label = paste("italic(b)[PLB]==",(round(PLB.bMLE.site.b[i],2))), parse = T) +
  annotate("text", x = 10, y = 3, label = paste("n =" ,(length(sitedata$Weight)))) +
    theme_classic() + 
    theme(axis.title = element_blank())
  AICdf$Number[i] <- length(sitedata$Weight)
  AICdf$lllognorm[i] <- lnormAIC(x)$lllognorm
  AICdf$AIClognorm[i] <- lnormAIC(x)$AIClognorm
  AICdf$llBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.biomass, xmax = siteinput$max.biomass, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$llBPL
  AICdf$AICBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.biomass, xmax = siteinput$max.biomass, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$AICBPL
}

leftlabel <- grid::textGrob(expression(paste("Number of birds with sizes", " ">=" ", italic("x"), "    ")), rot = 90)
bottomlabel <- grid::textGrob(expression(paste("Bird Weight, ", italic("x"), ~(g))))

ggsave(
  filename = "birdweights.pdf", 
  plot = marrangeGrob(grobs= siteb_plot, nrow=2, ncol=2,
                      left = leftlabel,
                      bottom = bottomlabel,
                      layout_matrix = rbind(c(1,2), c(3,4))), 
  width = 15, height = 9
)


write.csv(AICdf,'bird_AIC.csv')
write.csv(gof,'birdgof.csv')


# qqplots saved as pdf
pdf("qqplots_bird.pdf", width = 12, height = 6)
for(i in 1:nsites){
  par(mfrow=c(1,2))
  qqplot(FXinv(u = ppoints(siteinput$n), b = PLB.bMLE.site.b[i], xmin = min(sitex),
                             xmax = max(sitex)), x , xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(A)", sites[i], ": Bounded power law Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
  #qqline(x, distribution = function(p){
  #  FXinv(p, b = PLB.bMLE.site.b[i], xmin = min(sitex), xmax = max(sitex))
  #}, untf=T)
  qqplot(qlnorm(p = ppoints(siteinput$n), meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog), x, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(B)", sites[i], ": Log-normal Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
  #  qqline(x, distribution = function(p){qlnorm(p, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog)}, untf=T)
}
dev.off()

# # # saves the temp images from the plotting envrionment
# plots.dir.path <- list.files(tempdir(), pattern = "rs-graphics", full.names = TRUE)
# plots.png.paths <- list.files(plots.dir.path, pattern = ".png", full.names = TRUE)
# file.copy(from = plots.png.paths, to = "C:/Users/624225/OneDrive - hull.ac.uk/_BoxData/PhD/testingdistributions/qqplots/fish_logged")
