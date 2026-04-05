rm(list = ls())
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)

source("coralsizedistfuncs.R")

#load data
shark <- read.csv("SMART2324.csv")
shark <- shark %>% drop_na(Size..cm.)  # to remove NA
shark <- shark %>% drop_na(Bioregion)
shark <- shark %>% filter(!Species %in% c("Loggerhead Turtle", "Seal", "Tailor", "Snapper", "Black Marlin"))

sites <- unique(shark$Bioregion) # into 3 bioregions
nsites <- length(sites)
sites2 <- unique(shark$SD.region) # into 20 admin regions
nsites2 <- length(sites2)

# create empty lists to store things in - for shark lengths, 3 bioregions
PLB.bMLE.site.b <- numeric(nsites)
siteb_plot <- list()
bplqq <- list()
AICdf <- data.frame(Site = sites, Number = NA, llBPL = NA, lllognorm = NA, AICBPL = NA, AIClognorm = NA)
sigmadf <- data.frame(site = sites, lognorm = NA)
gof <- data.frame(site = sites, lognormX2 = NA, lognormdf = NA, lognormP = NA, BPLX2 = NA, BPLdf = NA, BPLP = NA)

# # run analysis for shark lengths by the three bioregions
# for(i in 1:nsites){
#   s <- sites[i]
#   sitedata <- shark %>% filter(Bioregion == s)
#   siteinput <- set.shark.params(sitedata$Size..cm.)
#   bML <- mle_b(Site == s, x = siteinput$biomass, sum_log_x = siteinput$sum.log.biomass,
#                x_min = siteinput$min.biomass, x_max = siteinput$max.biomass)
#   PLB.bMLE.site.b[i] <- bML[[1]] 
#   thetalnorm <- estimatelognormal(x = sitedata$Size..cm.)
#   x <- sitedata$Size..cm.
#   sitex = seq(min(sitedata$Size..cm.), max(sitedata$Size..cm.), length = 10000)
#   par(mfrow=c(1,2))
#   bplqq[[i]] <- qqplot(FXinv(u = ppoints(siteinput$n), b = PLB.bMLE.site.b[i], xmin = min(sitex),
#                              xmax = max(sitex)), x , xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(A)", sites[i], ": Power law Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
#   #qqline(x, distribution = function(p){
#   #  FXinv(p, b = PLB.bMLE.site.b[i], xmin = min(sitex), xmax = max(sitex))
#   #}, untf=T)
#   qqplot(qlnorm(p = ppoints(siteinput$n), meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog), x, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(B)", sites[i], ": Log-normal Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
#   #  qqline(x, distribution = function(p){qlnorm(p, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog)}, untf=T)
#   #rank plot   
#   sitey.PLB <- (1 - pPLB(x = sitex, b = PLB.bMLE.site.b[i], xmin = min(sitex),
#                          xmax = max(sitex))) * length(sitedata$Size..cm.)
#   sitey.lnorm <- plnorm(q = sitex, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog, lower.tail = FALSE, log.p = FALSE) * length(sitedata$Size..cm.)
#   siteb_plot[[i]] <- ggplot() +
#     geom_point(aes_(x = (sort(sitedata$Size..cm., decreasing=TRUE)), y = (1:length(sitedata$Size..cm.))),
#                color = "#666666", size = 2, alpha = 0.3) +
#     scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000), 
#                        limits = c(0.25, max(table(shark$Bioregion)))) +
#     scale_x_continuous(trans = 'log10',#breaks = c(-10,0,1,10,100),
#                        limits = c(50,max(sitedata$Size..cm.)))+
#     geom_line(aes_(x = sitex, y = sitey.PLB), col = 'black', lwd = 1) +
#     geom_line(aes_(x = sitex, y = sitey.lnorm), col = '#1B9E77', lwd = 1) +
#     labs(tag = LETTERS[i]) +
#    annotate("text", x = 100, y = 20, label = s) +
#     annotate("text", x = 100, y = 10, label = paste("italic(b)[PLB]==",(round(PLB.bMLE.site.b[i],2))), parse = T) +
#     annotate("text", x = 100, y =5, label = paste("n =" ,(length(sitedata$Size..cm.)))) +
#     theme_classic() + 
#     theme(axis.title = element_blank())
#   AICdf$Number[i] <- length(sitedata$Size..cm.)
#   AICdf$lllognorm[i] <- lnormAIC(x)$lllognorm
#   AICdf$AIClognorm[i] <- lnormAIC(x)$AIClognorm
#   AICdf$llBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.biomass, xmax = siteinput$max.biomass, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$llBPL
#   AICdf$AICBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.biomass, xmax = siteinput$max.biomass, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$AICBPL
# }
# 
# leftlabel <- grid::textGrob(expression(paste("Number of sharks with sizes", " ">=" ", italic("x"), "    ")), rot = 90)
# bottomlabel <- grid::textGrob(expression(paste("Shark length, ", italic("x"), ~(cm))))
# 
#  ggsave(
#    filename = "bioregionlengths.pdf", 
#    plot = marrangeGrob(grobs= siteb_plot, nrow=2, ncol=2,
#                        left = leftlabel,
#                        bottom = bottomlabel, 
#                        width = 15, height = 9))
#  
#  
# # 20 locations
#  PLB.bMLE.site2.b <- numeric(nsites2)
#  site2b_plot <- list()
#  bplqq2 <- list()
#  AICdf2 <- data.frame(Site = sites2, Number = NA, llBPL = NA, lllognorm = NA, AICBPL = NA, AIClognorm = NA)
#  
# # run analysis for shark lengths by the twenty locations
#  for(i in 1:nsites2){
#    s2 <- sites2[i]
#    sitedata2 <- shark %>% filter(SD.region == s2)
#    siteinput2 <- set.shark.params(sitedata2$Size..cm.)
#    bML2 <- mle_b(Site == s2, x = siteinput2$biomass, sum_log_x = siteinput2$sum.log.biomass,
#                 x_min = siteinput2$min.biomass, x_max = siteinput2$max.biomass)
#    PLB.bMLE.site2.b[i] <- bML2[[1]] 
#    thetalnorm2 <- estimatelognormal(x = sitedata2$Size..cm.)
#    x2 <- sitedata2$Size..cm.
#    sitex2 = seq(min(sitedata2$Size..cm.), max(sitedata2$Size..cm.), length = 10000)
#    par(mfrow=c(1,2))
#    bplqq2[[i]] <- qqplot(FXinv(u = ppoints(siteinput2$n), b = PLB.bMLE.site2.b[i], xmin = min(sitex2),
#                               xmax = max(sitex2)), x2 , xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(A)", sites2[i], ": Power law Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
#    #qqline(x, distribution = function(p){
#    #  FXinv(p, b = PLB.bMLE.site.b[i], xmin = min(sitex), xmax = max(sitex))
#    #}, untf=T)
#    qqplot(qlnorm(p = ppoints(siteinput2$n), meanlog = thetalnorm2$meanlog, sdlog = thetalnorm2$sdlog), x, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(B)", sites2[i], ": Log-normal Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
#    #  qqline(x, distribution = function(p){qlnorm(p, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog)}, untf=T)
#    #rank plot   
#    site2y.PLB <- (1 - pPLB(x = sitex2, b = PLB.bMLE.site2.b[i], xmin = min(sitex2),
#                           xmax = max(sitex2))) * length(sitedata2$Size..cm.)
#    site2y.lnorm <- plnorm(q = sitex2, meanlog = thetalnorm2$meanlog, sdlog = thetalnorm2$sdlog, lower.tail = FALSE, log.p = FALSE) * length(sitedata2$Size..cm.)
#    site2b_plot[[i]] <- ggplot() +
#      geom_point(aes_(x = (sort(sitedata2$Size..cm., decreasing=TRUE)), y = (1:length(sitedata2$Size..cm.))),
#                 color = "#666666", size = 2, alpha = 0.3) +
#      scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000), 
#                         limits = c(0.25, max(table(shark$SD.region)))) +
#      scale_x_continuous(trans = 'log10',#breaks = c(-10,0,1,10,100),
#                         limits = c(50,max(sitedata2$Size..cm.)))+
#      geom_line(aes_(x = sitex2, y = site2y.PLB), col = 'black', lwd = 1) +
#      geom_line(aes_(x = sitex2, y = site2y.lnorm), col = '#1B9E77', lwd = 1) +
#      labs(tag = LETTERS[i]) +
#      annotate("text", x = 100, y = 10, size = 6, label = s2) +
#      annotate("text", x = 100, y = 5, size = 5, label = paste("italic(b)[PLB]==",(round(PLB.bMLE.site2.b[i],2))), parse = T) +
#      annotate("text", x = 100, y =1, size = 5, label = paste("n =" ,(length(sitedata2$Size..cm.)))) +
#      theme_classic() + 
#      theme(axis.title = element_blank())
#    AICdf2$Number[i] <- length(sitedata2$Size..cm.)
#    AICdf2$lllognorm[i] <- lnormAIC(x2)$lllognorm
#    AICdf2$AIClognorm[i] <- lnormAIC(x2)$AIClognorm
#    AICdf2$llBPL[i] <- BPLAIC(C = getC(xmin = siteinput2$min.biomass, xmax = siteinput2$max.biomass, b = PLB.bMLE.site2.b[i]), b = PLB.bMLE.site2.b[i], x = x2)$llBPL
#    AICdf2$AICBPL[i] <- BPLAIC(C = getC(xmin = siteinput2$min.biomass, xmax = siteinput2$max.biomass, b = PLB.bMLE.site2.b[i]), b = PLB.bMLE.site2.b[i], x = x2)$AICBPL
#  } 
#  
#  
# # regionlengths <- marrangeGrob(grobs= site2b_plot, nrow=2, ncol=2,
# #                layout_matrix = rbind(c(1,2), c(3,4)))
# #  
# 
# ggsave(
#   filename = "locationlengths.pdf", 
#   plot = marrangeGrob(grobs= site2b_plot, nrow=2, ncol=2,
#                       left = leftlabel,
#                       bottom = bottomlabel,
#                       layout_matrix = rbind(c(1,2), c(3,4))), 
#   width = 15, height = 9
# )
# 
# write.csv(AICdf,'shark_AIC.csv')
# write.csv(AICdf2,'shark_region_AIC.csv')

# now repeat for Biomass

# create empty lists to store things in 
PLB.bMLE.site.b_biomass <- numeric(nsites)
siteb_plot_biomass <- list()
bplqq_biomass <- list()
AICdf_biomass <- data.frame(site = sites, Number = NA, llBPL = NA, lllognorm = NA, AICBPL = NA, AIClognorm = NA)

for(i in 1:nsites){
  s <- sites[i]
  sitedata <- shark %>% filter(Bioregion == s)
  siteinput_biomass <- set.fish.params(sitedata$biomass_kg)
  bML_biomass <- mle_b(Site == s, x = siteinput_biomass$biomass, sum_log_x = siteinput_biomass$sum.log.biomass,
               x_min = siteinput_biomass$min.biomass, x_max = siteinput_biomass$max.biomass)
  PLB.bMLE.site.b_biomass[i] <- bML_biomass[[1]] 
  thetalnorm_biomass <- estimatelognormal(x = sitedata$biomass_kg)
  
  #goodness-of-fit test for lognormal
  lngof <- lnormgof(x = sitedata$biomass_kg, mu = thetalnorm_biomass$meanlog, sigma = thetalnorm_biomass$sdlog)
  gof$lognormX2[i] <- lngof$X2
  gof$lognormdf[i] <- lngof$df
  gof$lognormP[i] <- lngof$P
  
  #goodness-of-fit test for bounded power law
  bplgof <- BPLgof(x = sitedata$biomass_kg, b = bML_biomass[[1]], xmin = min(sitedata$biomass_kg), xmax = max(sitedata$biomass_kg))
  gof$BPLX2[i] <- bplgof$X2
  gof$BPLdf[i] <- bplgof$df
  gof$BPLP[i] <- bplgof$P
  
  x_biomass <- sitedata$biomass_kg
  sitex_biomass = seq(min(sitedata$biomass_kg), max(sitedata$biomass_kg), length = 10000)
  par(mfrow=c(1,2))
  bplqq_biomass[[i]] <- qqplot(FXinv(u = ppoints(siteinput_biomass$n), b = PLB.bMLE.site.b_biomass[i], xmin = min(sitex_biomass),
                             xmax = max(sitex_biomass)), x_biomass , xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(A)", sites[i], ": Power law Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
  #qqline(x, distribution = function(p){
  #  FXinv(p, b = PLB.bMLE.site.b_biomass[i], xmin = min(sitex), xmax = max(sitex))
  #}, untf=T)
  qqplot(qlnorm(p = ppoints(siteinput_biomass$n), meanlog = thetalnorm_biomass$meanlog, sdlog = thetalnorm_biomass$sdlog), x_biomass, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(B)", sites[i], ": Log-normal Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
  #  qqline(x, distribution = function(p){qlnorm(p, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog)}, untf=T)
  #rank plot   
  sitey.PLB_biomass <- (1 - pPLB(x = sitex_biomass, b = PLB.bMLE.site.b_biomass[i], xmin = min(sitex_biomass),
                         xmax = max(sitex_biomass))) * length(sitedata$biomass_kg)
  sitey.lnorm_biomass <- plnorm(q = sitex_biomass, meanlog = thetalnorm_biomass$meanlog, sdlog = thetalnorm_biomass$sdlog, lower.tail = FALSE, log.p = FALSE) * length(sitedata$biomass_kg)
  siteb_plot_biomass[[i]] <- ggplot() +
    geom_point(aes_(x = (sort(sitedata$biomass_kg, decreasing=TRUE)), y = (1:length(sitedata$biomass_kg))),
               color = "#666666", size = 2, alpha = 0.3) +
    scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000), 
                       limits = c(0.25, max(table(shark$Bioregion)))) +
    scale_x_continuous(trans = 'log10',#breaks = c(-10,0,1,10,100),
                       limits = c(1,max(sitedata$biomass_kg)))+
    geom_line(aes_(x = sitex_biomass, y = sitey.PLB_biomass), col = 'black', lwd = 1) +
    geom_line(aes_(x = sitex_biomass, y = sitey.lnorm_biomass), col = '#1B9E77', lwd = 1) +
    labs(tag = LETTERS[i]) +
    annotate("text", x = 100, y = 5, label = s) +
    annotate("text", x = 100, y = 2.5, label = paste("italic(b)[PLB]==",(round(PLB.bMLE.site.b_biomass[i],2))), parse = T) +
    annotate("text", x = 100, y = 1.5, label = paste("n =" ,(length(sitedata$biomass_kg)))) +
    theme_classic() + 
    theme(axis.title = element_blank())
  AICdf_biomass$Number[i] <- length(sitedata$biomass_kg)
  AICdf_biomass$lllognorm[i] <- lnormAIC(x_biomass)$lllognorm
  AICdf_biomass$AIClognorm[i] <- lnormAIC(x_biomass)$AIClognorm
  AICdf_biomass$llBPL[i] <- BPLAIC(C = getC(xmin = siteinput_biomass$min.biomass, xmax = siteinput_biomass$max.biomass, b = PLB.bMLE.site.b_biomass[i]), b = PLB.bMLE.site.b_biomass[i], x = x_biomass)$llBPL
  AICdf_biomass$AICBPL[i] <- BPLAIC(C = getC(xmin = siteinput_biomass$min.biomass, xmax = siteinput_biomass$max.biomass, b = PLB.bMLE.site.b_biomass[i]), b = PLB.bMLE.site.b_biomass[i], x = x_biomass)$AICBPL
}

leftlabel_biomass <- grid::textGrob(expression(paste("Number of sharks with sizes", " ">=" ", italic("x"), "    ")), rot = 90)
bottomlabel_biomass <- grid::textGrob(expression(paste("Shark biomass, ", italic("x"), ~(kg))))

ggsave(
  filename = "bioregionbiomass.pdf", 
  plot = marrangeGrob(grobs= siteb_plot_biomass, nrow=2, ncol=2,
                      left = leftlabel_biomass,
                      bottom = bottomlabel_biomass, 
                      width = 15, height = 9))

# # create empty lists to store things in 
# PLB.bMLE.site2.b_biomass <- numeric(nsites2)
# site2b_plot_biomass <- list()
# bplqq2_biomass <- list()
# AICdf2_biomass <- data.frame(site = sites2, Number = NA, llBPL = NA, lllognorm = NA, AICBPL = NA, AIClognorm = NA)
# 
# for(i in 1:nsites2){
#   s2 <- sites2[i]
#   sitedata2 <- shark %>% filter(SD.region == s2)
#   siteinput2_biomass <- set.fish.params(sitedata2$biomass_kg)
#   bML2_biomass <- mle_b(Site == s2, x = siteinput2_biomass$biomass, sum_log_x = siteinput2_biomass$sum.log.biomass,
#                        x_min = siteinput2_biomass$min.biomass, x_max = siteinput2_biomass$max.biomass)
#   PLB.bMLE.site2.b_biomass[i] <- bML2_biomass[[1]] 
#   thetalnorm_biomass2 <- estimatelognormal(x = sitedata2$biomass_kg)
#   x2_biomass <- sitedata2$biomass_kg
#   sitex2_biomass = seq(min(sitedata2$biomass_kg), max(sitedata2$biomass_kg), length = 10000)
#   par(mfrow=c(1,2))
#   bplqq2_biomass[[i]] <- qqplot(FXinv(u = ppoints(siteinput2_biomass$n), b = PLB.bMLE.site2.b_biomass[i], xmin = min(sitex2_biomass),
#                                      xmax = max(sitex2_biomass)), x2_biomass , xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(A)", sites2[i], ": Power law Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
#   #qqline(x, distribution = function(p){
#   #  FXinv(p, b = PLB.bMLE.site.b_biomass[i], xmin = min(sitex), xmax = max(sitex))
#   #}, untf=T)
#   qqplot(qlnorm(p = ppoints(siteinput2_biomass$n), meanlog = thetalnorm_biomass2$meanlog, sdlog = thetalnorm_biomass2$sdlog), x2_biomass, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(B)", sites2[i], ": Log-normal Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
#   #  qqline(x, distribution = function(p){qlnorm(p, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog)}, untf=T)
#   #rank plot   
#   site2y.PLB_biomass <- (1 - pPLB(x = sitex2_biomass, b = PLB.bMLE.site2.b_biomass[i], xmin = min(sitex2_biomass),
#                                  xmax = max(sitex2_biomass))) * length(sitedata2$biomass_kg)
#   site2y.lnorm_biomass <- plnorm(q = sitex2_biomass, meanlog = thetalnorm_biomass2$meanlog, sdlog = thetalnorm_biomass2$sdlog, lower.tail = FALSE, log.p = FALSE) * length(sitedata2$biomass_kg)
#   site2b_plot_biomass[[i]] <- ggplot() +
#     geom_point(aes_(x = (sort(sitedata2$biomass_kg, decreasing=TRUE)), y = (1:length(sitedata2$biomass_kg))),
#                color = "#666666", size = 2, alpha = 0.3) +
#     scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000), 
#                        limits = c(0.25, max(table(shark$SD.region)))) +
#     scale_x_continuous(trans = 'log10',#breaks = c(-10,0,1,10,100),
#                        limits = c(1,max(sitedata$biomass_kg)))+
#     geom_line(aes_(x = sitex2_biomass, y = site2y.PLB_biomass), col = 'black', lwd = 1) +
#     geom_line(aes_(x = sitex2_biomass, y = site2y.lnorm_biomass), col = '#1B9E77', lwd = 1) +
#     labs(tag = LETTERS[i]) +
#     annotate("text", x = 10, y = 5, label = s2) +
#     annotate("text", x = 10, y = 2.5, label = paste("italic(b)[PLB]==",(round(PLB.bMLE.site2.b_biomass[i],2))), parse = T) +
#     annotate("text", x = 10, y = 1.5, label = paste("n =" ,(length(sitedata2$biomass_kg)))) +
#     theme_classic() + 
#     theme(axis.title = element_blank())
#   AICdf2_biomass$Number[i] <- length(sitedata2$biomass_kg)
#   AICdf2_biomass$lllognorm[i] <- lnormAIC(x2_biomass)$lllognorm
#   AICdf2_biomass$AIClognorm[i] <- lnormAIC(x2_biomass)$AIClognorm
#   AICdf2_biomass$llBPL[i] <- BPLAIC(C = getC(xmin = siteinput2_biomass$min.biomass, xmax = siteinput2_biomass$max.biomass, b = PLB.bMLE.site2.b_biomass[i]), b = PLB.bMLE.site2.b_biomass[i], x =  x2_biomass)$llBPL
#   AICdf2_biomass$AICBPL[i] <- BPLAIC(C = getC(xmin = siteinput2_biomass$min.biomass, xmax = siteinput2_biomass$max.biomass, b = PLB.bMLE.site2.b_biomass[i]), b = PLB.bMLE.site2.b_biomass[i], x =  x2_biomass)$AICBPL
# }
# 
# leftlabel_biomass <- grid::textGrob(expression(paste("Number of sharks with sizes", " ">=" ", italic("x"), "    ")), rot = 90)
# bottomlabel_biomass <- grid::textGrob(expression(paste("Shark biomass, ", italic("x"), ~(kg))))
# 
# 
# ggsave(
#   filename = "locationbiomass.pdf", 
#   plot = marrangeGrob(grobs= site2b_plot_biomass, nrow=2, ncol=2,
#                       left = leftlabel_biomass,
#                       bottom = bottomlabel_biomass,
#                       layout_matrix = rbind(c(1,2), c(3,4))), 
#   width = 15, height = 9
# )


#write.csv(AICdf_biomass,'shark_biomass_AIC.csv')
write.csv(AICdf_biomass,'shark_biomass_region_AIC.csv')

# # # # saves the temp images from the plotting envrionment
# plots.dir.path <- list.files(tempdir(), pattern = "rs-graphics", full.names = TRUE)
# plots.png.paths <- list.files(plots.dir.path, pattern = ".png", full.names = TRUE)
# file.copy(from = plots.png.paths, to = "C:/Users/tvkx991/OneDrive - University of Leeds/PhD/testingdistributions/testingdistributions")
