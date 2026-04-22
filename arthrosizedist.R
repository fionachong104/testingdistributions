rm(list = ls())
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)

source("coralsizedistfuncs.R")

#load data
arthro<- read.csv("Kollross_sizes.csv")
#tree <- tree %>% drop_na(aboveGroundLiveBiomass_kilograms)  # to remove dead trees
# axisscores <- read.csv("axisscores.csv", row.names=1)
# axisscores <- axisscores[order(axisscores$PC1), ]
arthro <- arthro[!is.na(arthro$bodylength_mm),]
arthro$siteName = paste(arthro$treatment, arthro$species, arthro$year, sep = "_")
sites <- unique(arthro$siteName)
#sites <- sites[match(row.names(axisscores), sites)] #re-ordering based on increasing PC1 scores
nsites <- length(sites)


# create empty lists to store things in 
PLB.bMLE.site.b <- numeric(nsites)
siteb_plot <- list()
bplqq <- list()
AICdf <- data.frame(Treatment = sites, Number = NA, llBPL = NA, lllognorm = NA, AICBPL = NA, AIClognorm = NA)
gof <- data.frame(site = sites, lognormX2 = NA, lognormdf = NA, lognormP = NA, BPLX2 = NA, BPLdf = NA, BPLP = NA)

for(i in 1:nsites){
  s <- sites[i]
  sitedata <- arthro %>% filter(siteName == s)
  siteinput <- set.arthro.params(sitedata$bodylength_mm)
  bML <- mle_b(Site == s, x = siteinput$biomass, sum_log_x = siteinput$sum.log.biomass,
               x_min = siteinput$min.biomass, x_max = siteinput$max.biomass)
  PLB.bMLE.site.b[i] <- bML[[1]] 
  thetalnorm <- estimatelognormal(x = sitedata$bodylength_mm)
  
  #goodness-of-fit test for lognormal
  lngof <- lnormgof(x = sitedata$bodylength_mm, mu = thetalnorm$meanlog, sigma = thetalnorm$sdlog)
  gof$lognormX2[i] <- lngof$X2
  gof$lognormdf[i] <- lngof$df
  gof$lognormP[i] <- lngof$P
  
  #goodness-of-fit test for bounded power law
  bplgof <- BPLgof(x = sitedata$bodylength_mm, b = bML[[1]], xmin = min(sitedata$bodylength_mm), xmax = max(sitedata$bodylength_mm))
  gof$BPLX2[i] <- bplgof$X2
  gof$BPLdf[i] <- bplgof$df
  gof$BPLP[i] <- bplgof$P
  
  x <- sitedata$bodylength_mm
  sitex = seq(min(sitedata$bodylength_mm), max(sitedata$bodylength_mm), length = 10000)
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
                         xmax = max(sitex))) * length(sitedata$bodylength_mm)
  sitey.lnorm <- plnorm(q = sitex, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog, lower.tail = FALSE, log.p = FALSE) * length(sitedata$bodylength_mm)
  siteb_plot[[i]] <- ggplot() +
    geom_point(aes_(x = (sort(sitedata$bodylength_mm, decreasing=TRUE)), y = (1:length(sitedata$bodylength_mm))),
               color = "#666666", size = 2, alpha = 0.3) +
    scale_y_continuous(trans = 'log10', breaks = c(1,10,100,350), 
                       limits = c(0.25, max(table(arthro$siteName)))) +
    scale_x_continuous(trans = 'log10',#breaks = c(-10,0,1,10,100),
                       limits = range(arthro$bodylength_mm))+
    geom_line(aes_(x = sitex, y = sitey.PLB), col = 'black', lwd = 1) +
    geom_line(aes_(x = sitex, y = sitey.lnorm), col = '#1B9E77', lwd = 1) +
    labs(tag = paste0("G", i))+
     annotate("text", x = 1, y = 10, label = s) +
   annotate("text", x = 1, y = 5, label = paste("italic(b)[PLB]==",(round(PLB.bMLE.site.b[i],2))), parse = T) +
   annotate("text", x = 1, y =1, label = paste("n =" ,(length(sitedata$bodylength_mm)))) +
    theme_classic() + 
    theme(axis.title = element_blank())
  AICdf$Number[i] <- length(sitedata$bodylength_mm)
  AICdf$lllognorm[i] <- lnormAIC(x)$lllognorm
  AICdf$AIClognorm[i] <- lnormAIC(x)$AIClognorm
  AICdf$llBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.biomass, xmax = siteinput$max.biomass, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$llBPL
  AICdf$AICBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.biomass, xmax = siteinput$max.biomass, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$AICBPL
}

leftlabel <- grid::textGrob(expression(paste("Number of arthropods with sizes", " ">=" ", italic("x"), "    ")), rot = 90)
bottomlabel <- grid::textGrob(expression(paste("Arthropods body length, ", italic("x"), ~(mm))))

ggsave(
  filename = "arthrobodylength.pdf", 
  plot = marrangeGrob(grobs= siteb_plot, nrow=2, ncol=2,
                      left = leftlabel,
                      bottom = bottomlabel,
                      layout_matrix = rbind(c(1,2), c(3,4))), 
  width = 15, height = 9
)

write.csv(AICdf,'arthro_AIC.csv')
write.csv(gof, 'gofarthro.csv')


hist(gof$BPLP, main = "(A) GOF test bounded power law p-values", breaks = seq(min(gof$BPLP), max(gof$BPLP) + 0.05, by = 0.05),#, col = badfit,
     xlim = c(0,1), ylim = c(0,16))
hist(gof$lognormP, main = "(B) GOF test log-normal p-values", breaks = seq(min(gof$BPLP), max(gof$BPLP) + 0.05, by = 0.05),#, col = badfit,
     xlim = c(0,1), ylim = c(0,16))


# # # saves the temp images from the plotting envrionment
# plots.dir.path <- list.files(tempdir(), pattern = "rs-graphics", full.names = TRUE)
# plots.png.paths <- list.files(plots.dir.path, pattern = ".png", full.names = TRUE)
# file.copy(from = plots.png.paths, to = "C:/Users/624225/OneDrive - hull.ac.uk/_BoxData/PhD/testingdistributions/qqplots/fish_logged")
