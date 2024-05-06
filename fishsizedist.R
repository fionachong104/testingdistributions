rm(list = ls())
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)

source("coralsizedistfuncs.R")

library(readxl)
#fish_size_spectra_data <- read_excel("C:/Users/624225/OneDrive - hull.ac.uk/_BoxData/PhD/testingdistributions/Carvalho21_code/fish_size_spectra_data.xlsx")
fishbiomass <- read_delim("fishbiomass_oneyear.csv")
axisscores <- read.csv("axisscores.csv", row.names=1)
axisscores <- axisscores[order(axisscores$PC1), ]

#reformatting Maria's fish data
#starting with BR2010
# BR_2010 <- as_tibble(BR_2010)
# TN <- BR_2010 %>%
#   group_by(Fish) %>%
#   summarise(total = sum(Number)) # this calculates how many fish there are, by species
fishbiomass <- as_tibble(fishbiomass)
Check_fishbiomass <- fishbiomass %>%
  group_by(Fish) %>%
  summarise(total = sum(Number))
#create data frame with each individual fish taking a row.
##make note of whether it is inividual biomass or not - each biomass value should be of a single fish not with them all grouped
# BRdf <- BR_2010 %>% 
#   uncount(Number) %>% 
#   group_by(Fish, Size) 

allfish <- fishbiomass %>%
  uncount(Number) %>%
  group_by(Fish, Size)

allfish$individual_biomass_kg <- as.numeric(allfish$individual_biomass_kg)
allfish %>% mutate(individual_biomass_kg = ifelse(is.na(individual_biomass_kg), 0, individual_biomass_kg))


sites <- unique(allfish$Site) #carn vs herb
sites <- sites[match(row.names(axisscores), sites)] #re-ordering based on increasing PC1 scores
nsites <- length(sites)

# create empty lists to store things in 
PLB.bMLE.site.b <- numeric(nsites)
siteb_plot <- list()
bplqq <- list()
AICdf <- data.frame(site = sites, llBPL = NA, lllognorm = NA, AICBPL = NA, AIClognorm = NA)


for(i in 1:nsites){
  s <- sites[i]
  sitedata <- allfish %>% filter(Site == s)
  siteinput <- set.fish.params(sitedata$individual_biomass_kg)
  bML <- mle_b(Site == s, x = siteinput$biomass, log_x = sitedata$log.biomass, sum_log_x = siteinput$sum.log.biomass,
               x_min = siteinput$min.biomass, x_max = siteinput$max.biomass)
  PLB.bMLE.site.b[i] <- bML[[1]] 
  thetalnorm <- estimatelognormal(x = sitedata$individual_biomass_kg)
  x <- sitedata$individual_biomass_kg
  sitex = seq(min(sitedata$individual_biomass_kg), max(sitedata$individual_biomass_kg), length = 1000)
  par(mfrow=c(1,2))
  bplqq[[i]] <- qqplot(FXinv(u = ppoints(siteinput$n), b = PLB.bMLE.site.b[i], xmin = min(sitex),
                             xmax = max(sitex)), x , xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", main = paste("(A)", sites[i], ": Power law Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
  qqline(x, distribution = function(p){
    FXinv(p, b = PLB.bMLE.site.b[i], xmin = min(sitex), xmax = max(sitex))
  })
  qqplot(qlnorm(p = ppoints(siteinput$n), meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog), x, xlab = "theoretical quantiles", ylab = "Sample quantiles", main = paste("(B)", sites[i], ": Log-normal Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
  qqline(x, distribution = function(p){qlnorm(p, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog)})
  #rank plot   
  sitey.PLB <- (1 - pPLB(x = sitex, b = PLB.bMLE.site.b[i], xmin = min(sitex),
                         xmax = max(sitex))) * length(sitedata$individual_biomass_kg)
  sitey.lnorm <- plnorm(q = sitex, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog, lower.tail = FALSE, log.p = FALSE) * length(sitedata$individual_biomass_kg)
  siteb_plot[[i]] <- ggplot() +
    geom_point(aes_(x = (sort(sitedata$individual_biomass_kg, decreasing=TRUE)), y = (1:length(sitedata$individual_biomass_kg))),
               color = "#666666", size = 2, alpha = 0.3) +
    scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000), 
                       limits = c(0.25, max(table(allfish$Site)))) +
    scale_x_continuous( trans = 'log10',breaks = c(1,10,100),
                        limits = range(allfish$individual_biomass_kg))+
    geom_line(aes_(x = sitex, y = sitey.PLB), col = 'black', lwd = 1) +
    geom_line(aes_(x = sitex, y = sitey.lnorm), col = '#1B9E77', lwd = 1) +
    labs(tag = LETTERS[i]) +
    annotate("text", x = 0.02, y = 10, label = s) +
    annotate("text", x = 0.02, y = 3, label = bquote(paste(italic(b)[PLB]==.(round(PLB.bMLE.site.b[i],2))))) +
    annotate("text", x = 100, y =1000, label = bquote(n == .(length(sitedata$individual_biomass_kg)))) +
    theme_classic() + 
    theme(axis.title = element_blank())
  AICdf$lllognorm[i] <- lnormAIC(x)$lllognorm
  AICdf$AIClognorm[i] <- lnormAIC(x)$AIClognorm
  AICdf$llBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.biomass, xmax = siteinput$max.biomass, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$llBPL
  AICdf$AICBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.biomass, xmax = siteinput$max.biomass, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$AICBPL
}

leftlabel <- grid::textGrob(expression(paste("Number of fish with sizes", " ">=" ", italic("x"), "    ")), rot = 90)
#bottomlabel <- grid::textGrob(expression(paste("Fish biomass, ", italic("x"), ~(kg))))
bottomlabel <- grid::textGrob(expression(paste("Fish biomass, ", italic("x"), ~(kg))))

siteb_plot <- grid.arrange(grobs = siteb_plot, ncol = 4, 
                           left = leftlabel,
                           bottom = bottomlabel)



# #input Carvalho fish df
# fish.df <- fish_size_spectra_data
# # Remove observer "ch" (initials PS) from analysis due to different overall size spectrum slope with 
# # # other divers in Lombok. See methods section.
# fish.df <- fish.df %>% filter(observer != "ch")
# carn.df <- fish.df %>% filter(tp == "Carnivore")
# herb.df <- fish.df %>% filter(tp == "Herbivore")
# 
# 


# ##for carvalho
# 
# for(i in 1:nsites){
#   s <- sites[i]
#   sitedata <- fish.df %>% filter(region == s)
#   siteinput <- set.params(sitedata$size_cm)
#   bML <- mle_b(Site == s, x = siteinput$Area, log_x = sitedata$log.Area, sum_log_x = siteinput$sum.log.Area,
#                x_min = siteinput$min.Area, x_max = siteinput$max.Area)
#   PLB.bMLE.site.b[i] <- bML[[1]] 
#   thetalnorm <- estimatelognormal(x = sitedata$size_cm)
#   x <- sitedata$size_cm
#   sitex = seq(min(sitedata$size_cm), max(sitedata$size_cm), length = 1000)
#   par(mfrow=c(2,2))
#   bplqq[[i]] <- qqplot(FXinv(u = ppoints(siteinput$n), b = PLB.bMLE.site.b[i], xmin = min(sitex),
#                              xmax = max(sitex)), x , xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", main = paste("(A)", sites[i], ": Power law Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
#   qqline(x, distribution = function(p){
#     FXinv(p, b = PLB.bMLE.site.b[i], xmin = min(sitex), xmax = max(sitex))
#   })
#   qqplot(qlnorm(p = ppoints(siteinput$n), meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog), x, xlab = "theoretical quantiles", ylab = "sample quantiles", main = paste("(B)", sites[i], ": Log-normal Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
#   qqline(x, distribution = function(p){qlnorm(p, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog)})
#   #rank plot   
#   sitey.PLB <- (1 - pPLB(x = sitex, b = PLB.bMLE.site.b[i], xmin = min(sitex),
#                          xmax = max(sitex))) * length(sitedata$size_cm)
#   sitey.lnorm <- plnorm(q = sitex, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog, lower.tail = FALSE, log.p = FALSE) * length(sitedata$size_cm)
#   siteb_plot[[i]] <- ggplot() +
#     geom_point(aes_(x = (sort(sitedata$size_cm, decreasing=TRUE)), y = (1:length(sitedata$size_cm))),
#                color = "cadetblue", size = 2, alpha = 0.3) +
#     scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000), 
#                        limits = c(0.25, max(table(fish.df$region)))) +
#     scale_x_continuous( breaks = c(10,20, 30, 40, 50, 60, 70),
#                        limits = range(fish.df$size_cm))+
#     geom_line(aes_(x = sitex, y = sitey.PLB), col = 'black', lwd = 1) +
#     geom_line(aes_(x = sitex, y = sitey.lnorm), col = 'green', lwd = 1) +
#     labs(tag = LETTERS[i]) +
#     annotate("text", x = 1, y = 10, label = s) +
#     annotate("text", x = 1, y = 3, label = bquote(paste(italic(b)[PLB]==.(round(PLB.bMLE.site.b[i],2))))) +
#     annotate("text", x = 3, y = 500, label = bquote(n == .(length(sitedata$size_cm)))) +
#     theme_classic() + 
#     theme(axis.title = element_blank())
#   AICdf$lllognorm[i] <- lnormAIC(x)$lllognorm
#   AICdf$AIClognorm[i] <- lnormAIC(x)$AIClognorm
#   AICdf$llBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.Area, xmax = siteinput$max.Area, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$llBPL
#   AICdf$AICBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.Area, xmax = siteinput$max.Area, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$AICBPL
# }
# 
# leftlabel <- grid::textGrob(expression(paste("Number of fish with sizes", " ">=" ", italic("x"), "    ")), rot = 90)
# #bottomlabel <- grid::textGrob(expression(paste("Fish biomass, ", italic("x"), ~(kg))))
# bottomlabel <- grid::textGrob(expression(paste("Fish length, ", italic("x"), ~(cm))))
# 
# siteb_plot <- grid.arrange(grobs = siteb_plot, ncol = 3, 
#                            left = leftlabel,
#                            bottom = bottomlabel)

# saves the temp images from the plotting envrionment
plots.dir.path <- list.files(tempdir(), pattern = "rs-graphics", full.names = TRUE)
plots.png.paths <- list.files(plots.dir.path, pattern = ".png", full.names = TRUE)
file.copy(from = plots.png.paths, to = "C:/Users/624225/OneDrive - hull.ac.uk/_BoxData/PhD/testingdistributions")
