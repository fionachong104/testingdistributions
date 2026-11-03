rm(list = ls())
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(moments)
library(cowplot)

source("coralsizedistfuncs.R")

fishbiomass <- read.delim("fishbiomass_oneyear.csv")
axisscores <- read.csv("axisscores.csv", row.names=1)
axisscores <- axisscores[order(axisscores$PC1), ]

#reformatting Maria's fish data
fishbiomass <- as_tibble(fishbiomass)
Check_fishbiomass <- fishbiomass %>%
  group_by(Fish) %>%
  summarise(total = sum(Number))

allfish <- fishbiomass %>%
  uncount(Number) %>%
  group_by(Fish, Size)

allfish$individual_biomass_kg <- as.numeric(allfish$individual_biomass_kg)
allfish %>% mutate(individual_biomass_kg = ifelse(is.na(individual_biomass_kg), 0, individual_biomass_kg))

#extract summary statistics 
all_fish <- allfish%>%
  group_by(Site) %>%
  summarise(
    fish_count = length(Fish),
    mean_indiv_length = round(mean(Size), digits = 3),
    mean_indiv_biomass = round(mean(individual_biomass_kg), digits = 3),
    median_indiv_biomass  = round(median(individual_biomass_kg), digits = 3),
    mean_log_biomass = round(mean(log(individual_biomass_kg)), digits = 3),
    median_log_biomass  = round(median(log(individual_biomass_kg)), digits = 3),
    sd = round(sd(log(individual_biomass_kg)), digits = 3), 
    skewness = round(skewness(log(individual_biomass_kg)), digits = 3),
    kurtosis = round(kurtosis(log(individual_biomass_kg)), digits = 3))

sites <- unique(allfish$Site) #missing three sites in the siteorder spreasheet
sites <- sites[match(row.names(axisscores), sites)]  #re-ordering based on order on hab het data
nsites <- length(sites)

# create empty lists to store things in 
PLB.bMLE.site.b <- numeric(nsites)
PLB.minLL.site.b <- numeric(nsites)
lnormparams <- data.frame(meanlog = numeric(nsites), sdlog = numeric(nsites))
siteb_plot <- list()
AICdf <- data.frame(site = sites, Number = NA, llBPL = NA, lllognorm = NA, AICBPL = NA, AIClognorm = NA)
gof <- data.frame(site = sites, lognormX2 = NA, lognormdf = NA, lognormP = NA, BPLX2 = NA, BPLdf = NA, BPLP = NA)
confintervals <- data.frame(site = sites, site_b = NA, lwr_CI95 = NA, upr_CI95 = NA, x_min = NA, x_max = NA, meanlog = NA, meanloglower95 = NA, meanlogupper95 = NA, sdlog = NA, sdloglower95 = NA, sdlogupper95 = NA)


for(i in 1:nsites){
  s <- sites[i]
  sitedata <- allfish %>% filter(Site == s)
  siteinput <- set.fish.params(sitedata$individual_biomass_kg)
  bML <- mle_b(Site == s, x = siteinput$biomass, sum_log_x = siteinput$sum.log.biomass,
               x_min = siteinput$min.biomass, x_max = siteinput$max.biomass)
  PLB.bMLE.site.b[i] <- bML[[1]] 
  PLB.minLL.site.b[i] <- bML[[2]]$minimum
  thetalnorm <- estimatelognormal(x = sitedata$individual_biomass_kg)
  lnormparams$meanlog[i] <- thetalnorm$meanlog
  lnormparams$sdlog[i] <- thetalnorm$sdlog
  #goodness-of-fit test for lognormal
  lngof <- lnormgof(x = sitedata$individual_biomass_kg, mu = thetalnorm$meanlog, sigma = thetalnorm$sdlog)
  gof$lognormX2[i] <- round(lngof$X2,3)
  gof$lognormdf[i] <- lngof$df
  gof$lognormP[i] <- signif(lngof$P,3)
  #goodness-of-fit test for bounded power law
  bplgof <- BPLgof(x = sitedata$individual_biomass_kg, b = bML[[1]], xmin = min(sitedata$individual_biomass_kg), xmax = max(sitedata$individual_biomass_kg))
  gof$BPLX2[i] <- round(bplgof$X2,3)
  gof$BPLdf[i] <- bplgof$df
  gof$BPLP[i] <- signif(bplgof$P,3)
  
  x <- sitedata$individual_biomass_kg
  sitex = seq(min(sitedata$individual_biomass_kg), max(sitedata$individual_biomass_kg), length = 10000)
 
  #rank plot   
  sitey.PLB <- (1 - pPLB(x = sitex, b = PLB.bMLE.site.b[i], xmin = min(sitex),
                         xmax = max(sitex))) * length(sitedata$individual_biomass_kg)
  sitey.lnorm <- plnorm(q = sitex, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog, lower.tail = FALSE, log.p = FALSE) * length(sitedata$individual_biomass_kg)
  
  model_df <- data.frame(
    sitex = sitex,
    value = c(sitey.PLB, sitey.lnorm),
    model = rep(c("Bounded Power Law", "Log-normal"), each = length(sitex))
  )
  
  siteb_plot[[i]] <- ggplot() +
    geom_point(aes_(x = (sort(sitedata$individual_biomass_kg, decreasing=TRUE)), y = (1:length(sitedata$individual_biomass_kg))),
               color = "#666666", size = 2, alpha = 0.3) +
    geom_line(
      data = model_df,
      aes(x = sitex, y = value, colour = model),
      linewidth = 1
    ) +
    scale_colour_manual(
      values = c("Bounded Power Law" = "black", "Log-normal" = "#1B9E77"),
      name = "Model"
    ) +
    scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000), 
                       limits = c(0.25, max(table(allfish$Site)))) +
    scale_x_continuous(trans = 'log10', breaks = c(0.0001,0.01,1,100),
                        limits = range(allfish$individual_biomass_kg))+
    labs(tag = LETTERS[i]) +
    annotate("text", x = 0.002, y = 10, label = s) +
    annotate("text", x = 0.002, y = 5, label = paste("n =" ,(length(sitedata$individual_biomass_kg)))) +
    theme_classic() + 
    theme(axis.title = element_blank(), legend.position = "none")
  AICdf$Number[i] <- length(sitedata$individual_biomass_kg)
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

write.csv(AICdf,'fish_AIC.csv')
write.csv(gof, "fishgof.csv")
write.csv(confintervals, "fishparam.csv")

# getting the legend from the first siteb_plot
legend_plot <- get_legend(
  siteb_plot[[1]] +
    theme(legend.position = "top",
          legend.direction = "horizontal",
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16))
)

# putting the plots into a 5x4 grid

leftlabel <- grid::textGrob(expression(paste("Number of fish with sizes", " ">=" ", italic("x"), "    ")), rot = 90,
                            gp = grid::gpar(fontsize = 16))
bottomlabel <- grid::textGrob(expression(paste("Fish biomass, ", italic("x"), ~(kg))), gp = grid::gpar(fontsize = 16))

panel_grid <- arrangeGrob(
  grobs = siteb_plot,
  nrow  = 5,
  ncol  = 4,
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
  filename = "fishbiomass.pdf", 
  final_plot, 
  width = 15, height = 13
)

# qqplots saved as pdf
pdf("qqplots_fish.pdf", width = 12, height = 6)
for(i in 1:nsites){
  s <- sites[i]
  sitedata <- allfish %>% filter(Site == s)
  siteinput <- set.fish.params(sitedata$individual_biomass_kg)
  x <- sitedata$individual_biomass_kg
  sitex = seq(min(sitedata$individual_biomass_kg), max(sitedata$individual_biomass_kg), length = 10000)
  par(mfrow=c(1,2))
  qqplot(FXinv(u = ppoints(siteinput$n), b = PLB.bMLE.site.b[i], xmin = min(sitex),
                             xmax = max(sitex)), x , xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(A)", sites[i], ": Bounded Power law Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
  qqplot(qlnorm(p = ppoints(siteinput$n), meanlog = lnormparams$meanlog[i], sdlog = lnormparams$sdlog[i]), x, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(B)", sites[i], ": Log-normal"), pch = 16, col = adjustcolor("black", 0.25))
}
dev.off()
