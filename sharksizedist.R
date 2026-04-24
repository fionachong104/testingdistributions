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

# create empty lists to store things in 
PLB.bMLE.site.b_biomass <- numeric(nsites)
siteb_plot_biomass <- list()
AICdf_biomass <- data.frame(site = sites, Number = NA, llBPL = NA, lllognorm = NA, AICBPL = NA, AIClognorm = NA)
gof <- data.frame(site = sites, lognormX2 = NA, lognormdf = NA, lognormP = NA, BPLX2 = NA, BPLdf = NA, BPLP = NA)

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

leftlabel <- grid::textGrob(expression(paste("Number of sharks with sizes", " ">=" ", italic("x"), "    ")), rot = 90)
bottomlabel<- grid::textGrob(expression(paste("Shark biomass, ", italic("x"), ~(kg))))

ggsave(
  filename = "sharkbiomass.pdf", 
  plot = marrangeGrob(grobs= siteb_plot_biomass, nrow=2, ncol=2,
                      left = leftlabel,
                      bottom = bottomlabel,
                      layout_matrix = rbind(c(1,2), c(3,4))), 
  width = 15, height = 9
)

# qqplots saved as pdf
pdf("qqplots_sharks.pdf", width = 12, height = 6)
for(i in 1:nsites){
par(mfrow=c(1,2))
qqplot(FXinv(u = ppoints(siteinput_biomass$n), b = PLB.bMLE.site.b_biomass[i], xmin = min(sitex_biomass),
             xmax = max(sitex_biomass)), x_biomass , xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(A)", sites[i], ": Bounded power law Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
#qqline(x, distribution = function(p){
#  FXinv(p, b = PLB.bMLE.site.b_biomass[i], xmin = min(sitex), xmax = max(sitex))
#}, untf=T)
qqplot(qlnorm(p = ppoints(siteinput_biomass$n), meanlog = thetalnorm_biomass$meanlog, sdlog = thetalnorm_biomass$sdlog), x_biomass, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(B)", sites[i], ": Log-normal Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
#  qqline(x, distribution = function(p){qlnorm(p, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog)}, untf=T)
}
dev.off()

write.csv(AICdf_biomass,'shark_biomass_AIC.csv')
write.csv(gof,'sharkgof.csv')
