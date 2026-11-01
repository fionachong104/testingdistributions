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
PLB.minLL.site.b <- numeric(nsites)
siteb_plot <- list()
AICdf <- data.frame(site = sites, llBPL = NA, lllognorm = NA, AICBPL = NA, AIClognorm = NA)
gof <- data.frame(site = sites, lognormX2 = NA, lognormdf = NA, lognormP = NA, BPLX2 = NA, BPLdf = NA, BPLP = NA)
confintervals <- data.frame(site = sites, site_b = NA, lwr_CI95 = NA, upr_CI95 = NA, x_min = NA, x_max = NA, meanlog = NA, meanloglower95 = NA, meanlogupper95 = NA, sdlog = NA, sdloglower95 = NA, sdlogupper95 = NA)

for(i in 1:nsites){
  s <- sites[i]
  sitedata <- fish.df %>% filter(region == s)
  siteinput <- set.fish.params(sitedata$biomass_kg)
  bML <- mle_b(Site == s, x = siteinput$biomass, sum_log_x = siteinput$sum.log.biomass,
               x_min = siteinput$min.biomass, x_max = siteinput$max.biomass)
  PLB.bMLE.site.b[i] <- bML[[1]]
  PLB.minLL.site.b[i] <- bML[[2]]
  thetalnorm <- estimatelognormal(x = sitedata$biomass_kg)
  
  #goodness-of-fit test for lognormal
  lngof <- lnormgof(x = sitedata$biomass_kg, mu = thetalnorm$meanlog, sigma = thetalnorm$sdlog)
  gof$lognormX2[i] <- lngof$X2
  gof$lognormdf[i] <- lngof$df
  gof$lognormP[i] <- signif(lngof$P,3)
  
  #goodness-of-fit test for bounded power law
  bplgof <- BPLgof(x = sitedata$biomass_kg, b = bML[[1]], xmin = min(sitedata$biomass_kg), xmax = max(sitedata$biomass_kg))
  gof$BPLX2[i] <- bplgof$X2
  gof$BPLdf[i] <- bplgof$df
  gof$BPLP[i] <- signif(bplgof$P,3)
  
  x <- sitedata$biomass_kg
  sitex = seq(min(sitedata$biomass_kg), max(sitedata$biomass_kg), length = 1000)
 
  #rank plot
  sitey.PLB <- (1 - pPLB(x = sitex, b = PLB.bMLE.site.b[i], xmin = min(sitex),
                         xmax = max(sitex))) * length(sitedata$biomass_kg)
  sitey.lnorm <- plnorm(q = sitex, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog, lower.tail = FALSE, log.p = FALSE) * length(sitedata$biomass_kg)
  
  model_df <- data.frame(
    sitex = sitex,
    value = c(sitey.PLB, sitey.lnorm),
    model = rep(c("Bounded Power Law", "Log-normal"), each = length(sitex))) 
  
  siteb_plot[[i]] <- ggplot() +
    geom_line(
      data = model_df,
      aes(x = sitex, y = value, colour = model),
      linewidth = 1
    ) +
    scale_colour_manual(
      values = c("Bounded Power Law" = "black", "Log-normal" = "#1B9E77"),
      name = "Model"
    ) +
    geom_point(aes_(x = (sort(sitedata$biomass_kg, decreasing=TRUE)), y = (1:length(sitedata$biomass_kg))),
               color = "#666666", size = 2, alpha = 0.3) +
    scale_y_continuous(trans = 'log10', breaks = c(1,10,100,1000,10000),
                       limits = c(0.25, max(table(fish.df$region)))) +
    scale_x_continuous(trans = 'log10', breaks = c(0.01,0.5,1,5),
                        limits = range(fish.df$biomass_kg))+
    geom_line(aes_(x = sitex, y = sitey.PLB), col = 'black', lwd = 1) +
    geom_line(aes_(x = sitex, y = sitey.lnorm), col = '#1B9E77', lwd = 1) +
    labs(tag = LETTERS[i]) +
    annotate("text", x = 0.5, y = 10, label = s) +
    annotate("text", x = 0.5, y = 3, label = paste("n =" ,(length(sitedata$biomass_kg)))) +
    theme_classic() +
    theme(axis.title = element_blank(), legend.position = "none")
  AICdf$lllognorm[i] <- lnormAIC(x)$lllognorm
  AICdf$AIClognorm[i] <- lnormAIC(x)$AIClognorm
  AICdf$llBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.biomass, xmax = siteinput$max.biomass, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$llBPL
  AICdf$AICBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.biomass, xmax = siteinput$max.biomass, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$AICBPL
  confintervals$site_b[i] <- round(PLB.bMLE.site.b[i],3)
  confintervals$lwr_CI95[i] <- round(slope.conf.int(PLB.bMLE.site.b[[i]], PLB.minLL.site.b[[i]], siteinput)[1],3)
  confintervals$upr_CI95[i] <- round(slope.conf.int(PLB.bMLE.site.b[[i]], PLB.minLL.site.b[[i]], siteinput)[2],3)
  confintervals$x_min[i] <- round(siteinput$min.biomass, 5)
  confintervals$x_max[i] <- round(siteinput$max.biomass, 3)
  confintervals$meanlog[i] <- round(thetalnorm$meanlog, 3)
  confintervals$meanloglower95[i] <- round(thetalnorm$meanlower, 3)
  confintervals$meanlogupper95[i] <- round(thetalnorm$meanupper, 3)
  confintervals$sdlog[i] <- round(thetalnorm$sdlog, 3)
  confintervals$sdloglower95[i] <- round(thetalnorm$sdlower, 3)
  confintervals$sdlogupper95[i] <- round(thetalnorm$sdupper, 3)
  
}

leftlabel <- grid::textGrob(expression(paste("Number of fish with sizes", " ">=" ", italic("x"), "    ")), rot = 90)
bottomlabel <- grid::textGrob(expression(paste("Fish biomass, ", italic("x"), ~(kg))))
#bottomlabel <- grid::textGrob(expression(paste("Fish length, ", italic("x"), ~(cm))))

# getting the legend from the first siteb_plot
legend_plot <- get_legend(
  siteb_plot[[1]] +
    theme(legend.position = "top",
          legend.direction = "horizontal",
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16))
)

panel_grid <- arrangeGrob(
  grobs = siteb_plot,
  nrow  = 2,
  ncol  = 2,
  left  = leftlabel,
  bottom = bottomlabel
)

final_plot <- plot_grid(
  legend_plot,
  panel_grid,
  ncol = 1,
  rel_heights = c(0.05, 1)
)


ggsave(
  filename = "carvalhobiomass.pdf", 
  final_plot, 
  width = 15, height = 13
)

write.csv(AICdf,'carvalho_AIC.csv')
write.csv(gof,'carvalhogof.csv')
write.csv(confintervals, "carvalhoparam.csv")


# qqplots saved as pdf
pdf("qqplots_carvalho.pdf", width = 12, height = 6)
for(i in 1:nsites){
  par(mfrow=c(1,2))
  qqplot(FXinv(u = ppoints(siteinput$n), b = PLB.bMLE.site.b[i], xmin = min(sitex),
               xmax = max(sitex)), x , xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(A)", sites[i], ": Bounded power law"), pch = 16, col = adjustcolor("black", 0.25))
  #qqline(x, distribution = function(p){
  #   FXinv(p, b = PLB.bMLE.site.b[i], xmin = min(sitex), xmax = max(sitex))
  #}, untf=T)
  qqplot(qlnorm(p = ppoints(siteinput$n), meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog), x, xlab = "Theoretical quantiles", ylab = "Sample quantiles", log = "xy", main = paste("(B)", sites[i], ": Log-normal"), pch = 16, col = adjustcolor("black", 0.25))
  #  qqline(x, distribution = function(p){qlnorm(p, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog)}, untf=T)
  }
dev.off()
