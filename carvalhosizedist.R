rm(list = ls())
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(readxl)

HELLO

source("coralsizedistfuncs.R")

#input Carvalho fish df
fish_size_spectra_data <- read_excel("fish_size_spectra_data.xlsx")
fish.df <- fish_size_spectra_data
# Remove observer "ch" (initials PS) from analysis due to different overall size spectrum slope with
# # other divers in Lombok. See methods section.
fish.df <- fish.df %>% filter(observer != "ch")
carn.df <- fish.df %>% filter(tp == "Carnivore")
herb.df <- fish.df %>% filter(tp == "Herbivore")

sites <- unique(fish.df$region) #carn vs herb
nsites <- length(sites)

# create empty lists to store things in 
PLB.bMLE.site.b <- numeric(nsites)
siteb_plot <- list()
bplqq <- list()
AICdf <- data.frame(site = sites, llBPL = NA, lllognorm = NA, AICBPL = NA, AIClognorm = NA)



##for carvalho

for(i in 1:nsites){
  s <- sites[i]
  sitedata <- fish.df %>% filter(region == s)
  siteinput <- set.fish.params(sitedata$biomass_kg)
  bML <- mle_b(Site == s, x = siteinput$biomass, log_x = sitedata$log.biomass, sum_log_x = siteinput$sum.log.biomass,
               x_min = siteinput$min.biomass, x_max = siteinput$max.biomass)
  PLB.bMLE.site.b[i] <- bML[[1]]
  thetalnorm <- estimatelognormal(x = sitedata$biomass_kg)
  x <- sitedata$biomass_kg
  sitex = seq(min(sitedata$biomass_kg), max(sitedata$biomass_kg), length = 1000)
  par(mfrow=c(1,2))
  bplqq[[i]] <- qqplot(FXinv(u = ppoints(siteinput$n), b = PLB.bMLE.site.b[i], xmin = min(sitex),
                             xmax = max(sitex)), x , xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", main = paste("(A)", sites[i], ": Power law Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
  qqline(x, distribution = function(p){
    FXinv(p, b = PLB.bMLE.site.b[i], xmin = min(sitex), xmax = max(sitex))
  })
  qqplot(qlnorm(p = ppoints(siteinput$n), meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog), x, xlab = "theoretical quantiles", ylab = "sample quantiles", main = paste("(B)", sites[i], ": Log-normal Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
  qqline(x, distribution = function(p){qlnorm(p, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog)})
  #rank plot
  sitey.PLB <- (1 - pPLB(x = sitex, b = PLB.bMLE.site.b[i], xmin = min(sitex),
                         xmax = max(sitex))) * length(sitedata$biomass_kg)
  sitey.lnorm <- plnorm(q = sitex, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog, lower.tail = FALSE, log.p = FALSE) * length(sitedata$biomass_kg)
  siteb_plot[[i]] <- ggplot() +
    geom_point(aes_(x = (sort(sitedata$biomass_kg, decreasing=TRUE)), y = (1:length(sitedata$biomass_kg))),
               color = "cadetblue", size = 2, alpha = 0.3) +
    scale_y_continuous(trans = 'log10', breaks = c(1,10,100,1000,10000),
                       limits = c(0.25, max(table(fish.df$region)))) +
    scale_x_continuous( breaks = c(0,1,5,10),
                        limits = range(fish.df$biomass_kg))+
    geom_line(aes_(x = sitex, y = sitey.PLB), col = 'black', lwd = 1) +
    geom_line(aes_(x = sitex, y = sitey.lnorm), col = 'green', lwd = 1) +
    labs(tag = LETTERS[i]) +
    annotate("text", x = 1.5, y = 10, label = s) +
    annotate("text", x = 1.5, y = 3, label = bquote(paste(italic(b)[PLB]==.(round(PLB.bMLE.site.b[i],2))))) +
    annotate("text", x = 4, y = 1000, label = bquote(n == .(length(sitedata$biomass_kg)))) +
    theme_classic() +
    theme(axis.title = element_blank())
  AICdf$lllognorm[i] <- lnormAIC(x)$lllognorm
  AICdf$AIClognorm[i] <- lnormAIC(x)$AIClognorm
  AICdf$llBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.biomass, xmax = siteinput$max.biomass, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$llBPL
  AICdf$AICBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.biomass, xmax = siteinput$max.biomass, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$AICBPL
}

leftlabel <- grid::textGrob(expression(paste("Number of fish with sizes", " ">=" ", italic("x"), "    ")), rot = 90)
bottomlabel <- grid::textGrob(expression(paste("Fish biomass, ", italic("x"), ~(kg))))
#bottomlabel <- grid::textGrob(expression(paste("Fish length, ", italic("x"), ~(cm))))

siteb_plot <- grid.arrange(grobs = siteb_plot, ncol = 3,
                           left = leftlabel,
                           bottom = bottomlabel)


write.csv(AICdf,'carvalho_AIC.csv')
