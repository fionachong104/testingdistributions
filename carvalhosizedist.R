rm(list = ls())
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(readxl)

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
gof <- data.frame(site = sites, lognormX2 = NA, lognormdf = NA, lognormP = NA, BPLX2 = NA, BPLdf = NA, BPLP = NA)



##for carvalho

for(i in 1:nsites){
  s <- sites[i]
  sitedata <- fish.df %>% filter(region == s)
  siteinput <- set.fish.params(sitedata$biomass_kg)
  bML <- mle_b(Site == s, x = siteinput$biomass, sum_log_x = siteinput$sum.log.biomass,
               x_min = siteinput$min.biomass, x_max = siteinput$max.biomass)
  PLB.bMLE.site.b[i] <- bML[[1]]
  thetalnorm <- estimatelognormal(x = sitedata$biomass_kg)
  
  #goodness-of-fit test for lognormal
  lngof <- lnormgof(x = sitedata$biomass_kg, mu = thetalnorm$meanlog, sigma = thetalnorm$sdlog)
  gof$lognormX2[i] <- lngof$X2
  gof$lognormdf[i] <- lngof$df
  gof$lognormP[i] <- lngof$P
  
  #goodness-of-fit test for bounded power law
  bplgof <- BPLgof(x = sitedata$biomass_kg, b = bML[[1]], xmin = min(sitedata$biomass_kg), xmax = max(sitedata$biomass_kg))
  gof$BPLX2[i] <- bplgof$X2
  gof$BPLdf[i] <- bplgof$df
  gof$BPLP[i] <- bplgof$P
  
  x <- sitedata$biomass_kg
  sitex = seq(min(sitedata$biomass_kg), max(sitedata$biomass_kg), length = 1000)
  par(mfrow=c(1,2))
  bplqq[[i]] <- qqplot(FXinv(u = ppoints(siteinput$n), b = PLB.bMLE.site.b[i], xmin = min(sitex),
                             xmax = max(sitex)), x , xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(A)", sites[i], ": Power law"), pch = 16, col = adjustcolor("black", 0.25))
  #qqline(x, distribution = function(p){
 #   FXinv(p, b = PLB.bMLE.site.b[i], xmin = min(sitex), xmax = max(sitex))
  #}, untf=T)
  qqplot(qlnorm(p = ppoints(siteinput$n), meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog), x, xlab = "Theoretical quantiles", ylab = "Sample quantiles", log = "xy", main = paste("(B)", sites[i], ": Log-normal"), pch = 16, col = adjustcolor("black", 0.25))
#  qqline(x, distribution = function(p){qlnorm(p, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog)}, untf=T)
  #rank plot
  sitey.PLB <- (1 - pPLB(x = sitex, b = PLB.bMLE.site.b[i], xmin = min(sitex),
                         xmax = max(sitex))) * length(sitedata$biomass_kg)
  sitey.lnorm <- plnorm(q = sitex, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog, lower.tail = FALSE, log.p = FALSE) * length(sitedata$biomass_kg)
  siteb_plot[[i]] <- ggplot() +
    geom_point(aes_(x = (sort(sitedata$biomass_kg, decreasing=TRUE)), y = (1:length(sitedata$biomass_kg))),
               color = "#666666", size = 2, alpha = 0.3) +
    scale_y_continuous(trans = 'log10', breaks = c(1,10,100,1000,10000),
                       limits = c(0.25, max(table(fish.df$region)))) +
    scale_x_continuous( breaks = c(0.01,0.5,1,5),
                        limits = range(fish.df$biomass_kg))+
    geom_line(aes_(x = sitex, y = sitey.PLB), col = 'black', lwd = 1) +
    geom_line(aes_(x = sitex, y = sitey.lnorm), col = '#1B9E77', lwd = 1) +
    labs(tag = paste0("C", i)) +
    annotate("text", x = 1.5, y = 10, label = s) +
    annotate("text", x = 1.5, y = 5, label = paste("italic(b)[PLB]==",(round(PLB.bMLE.site.b[i],2))), parse = T)  +
    annotate("text", x = 1.5, y = 3, label = paste("n =" ,(length(sitedata$biomass_kg)))) +
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

ggsave(
  filename = "carvalhobiomass.pdf", 
  plot = marrangeGrob(grobs= siteb_plot, nrow=2, ncol=2,
                      left = leftlabel,
                      bottom = bottomlabel,
                      layout_matrix = rbind(c(1,2), c(3,4))), 
  width = 15, height = 9
)


write.csv(AICdf,'carvalho_AIC.csv')
write.csv(gof,'carvalhogof.csv')

hist(gof$BPLP, main = "(A) GOF test bounded power law p-values", breaks = seq(min(gof$BPLP), max(gof$BPLP) + 0.05, by = 0.05),#, col = badfit,
     xlim = c(0,1), ylim = c(0,20))
hist(gof$lognormP, main = "(B) GOF test log-normal p-values", breaks = seq(min(gof$lognormP), max(gof$lognormP) + 0.05, by = 0.05), #col = badfit,
     xlim = c(0,1), ylim = c(0,20))