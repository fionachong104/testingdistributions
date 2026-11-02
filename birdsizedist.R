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
PLB.minLL.site.b <- numeric(nsites)
siteb_plot <- list()
AICdf <- data.frame(Site = sites, Number = NA, llBPL = NA, lllognorm = NA, AICBPL = NA, AIClognorm = NA)
gof <- data.frame(site = sites, lognormX2 = NA, lognormdf = NA, lognormP = NA, BPLX2 = NA, BPLdf = NA, BPLP = NA)
confintervals <- data.frame(site = sites, site_b = NA, lwr_CI95 = NA, upr_CI95 = NA, x_min = NA, x_max = NA, meanlog = NA, meanloglower95 = NA, meanlogupper95 = NA, sdlog = NA, sdloglower95 = NA, sdlogupper95 = NA)


for(i in 1:nsites){
  s <- sites[i]
  sitedata <- bird %>% filter(siteName == s)
  siteinput <- set.bird.params(sitedata$Weight)
  bML <- mle_b(Site == s, x = siteinput$biomass, sum_log_x = siteinput$sum.log.biomass,
               x_min = siteinput$min.biomass, x_max = siteinput$max.biomass)
  PLB.bMLE.site.b[i] <- bML[[1]]
  PLB.minLL.site.b[i] <- bML[[2]]$minimum
  thetalnorm <- estimatelognormal(x = sitedata$Weight)

  #goodness-of-fit test for lognormal
  lngof <- lnormgof(x = sitedata$Weight, mu = thetalnorm$meanlog, sigma = thetalnorm$sdlog)
  gof$lognormX2[i] <- round(lngof$X2,3)
  gof$lognormdf[i] <- lngof$df
  gof$lognormP[i] <- signif(lngof$P,3)

  #goodness-of-fit test for bounded power law
  bplgof <- BPLgof(x = sitedata$Weight, b = bML[[1]], xmin = min(sitedata$Weight), xmax = max(sitedata$Weight))
  gof$BPLX2[i] <- round(bplgof$X2,3)
  gof$BPLdf[i] <- bplgof$df
  gof$BPLP[i] <- signif(bplgof$P,3)

  x <- sitedata$Weight
  sitex = seq(min(sitedata$Weight), max(sitedata$Weight), length = 10000)
    #rank plot   
  sitey.PLB <- (1 - pPLB(x = sitex, b = PLB.bMLE.site.b[i], xmin = min(sitex),
                         xmax = max(sitex))) * length(sitedata$Weight)
  sitey.lnorm <- plnorm(q = sitex, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog, lower.tail = FALSE, log.p = FALSE) * length(sitedata$Weight)
  
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
  annotate("text", x = 10, y = 5, label = paste("n =" ,(length(sitedata$Weight)))) +
    theme_classic() + 
    theme(axis.title = element_blank(), legend.position = "none")
  AICdf$Number[i] <- length(sitedata$Weight)
  AICdf$lllognorm[i] <- lnormAIC(x)$lllognorm
  AICdf$AIClognorm[i] <- lnormAIC(x)$AIClognorm
  AICdf$llBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.biomass, xmax = siteinput$max.biomass, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$llBPL
  AICdf$AICBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.biomass, xmax = siteinput$max.biomass, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$AICBPL
  confintervals$site_b[i] <- round(PLB.bMLE.site.b[i],3)
  confintervals$lwr_CI95[i] <- round(slope.conf.int(PLB.bMLE.site.b[[i]], PLB.minLL.site.b[[i]], siteinput)[1],3)
  confintervals$upr_CI95[i] <- round(slope.conf.int(PLB.bMLE.site.b[[i]], PLB.minLL.site.b[[i]], siteinput)[2],3)
  confintervals$x_min[i] <- round(siteinput$min.biomass, 3)
  confintervals$x_max[i] <- round(siteinput$max.biomass, 3)
  confintervals$meanlog[i] <- round(thetalnorm$meanlog, 3)
  confintervals$meanloglower95[i] <- round(thetalnorm$meanlower, 3)
  confintervals$meanlogupper95[i] <- round(thetalnorm$meanupper, 3)
  confintervals$sdlog[i] <- round(thetalnorm$sdlog, 3)
  confintervals$sdloglower95[i] <- round(thetalnorm$sdlower, 3)
  confintervals$sdlogupper95[i] <- round(thetalnorm$sdupper, 3)
  
}

leftlabel <- grid::textGrob(expression(paste("Number of birds with sizes", " ">=" ", italic("x"), "    ")), rot = 90)
bottomlabel <- grid::textGrob(expression(paste("Bird Weight, ", italic("x"), ~(g))))

# getting the legend from the first siteb_plot
legend_plot <- get_legend(
  siteb_plot[[1]] +
    theme(legend.position = "top",
          legend.direction = "horizontal",
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16))
)

# Split plots into groups of 4
plot_pages <- split(
  siteb_plot,
  ceiling(seq_along(siteb_plot) / 4)
)

pages <- lapply(plot_pages, function(x){
  
  panel_grid <- arrangeGrob(
    grobs = x,
    nrow = 2,
    ncol = 2,
    left = leftlabel,
    bottom = bottomlabel
  )
  
  plot_grid(
    legend_plot,
    panel_grid,
    ncol = 1,
    rel_heights = c(0.05, 1)
  )
})

finalplot <- marrangeGrob(
  pages,
  nrow = 1,
  ncol = 1,
  top = quote(paste("Page", g, "of", npages))
)

ggsave(
  "birdweights.pdf",
  finalplot,
  width = 12,
  height = 8.5,
  units = "in"
)


write.csv(AICdf,'bird_AIC.csv')
write.csv(gof,'birdgof.csv')
write.csv(confintervals, "birdparam.csv")

# qqplots saved as pdf
pdf("qqplots_bird.pdf", width = 12, height = 6)
for(i in 1:nsites){
  s <- sites[i]
  sitedata <- bird %>% filter(siteName == s)
  siteinput <- set.bird.params(sitedata$Weight)
  x <- sitedata$Weight
  sitex = seq(min(sitedata$Weight), max(sitedata$Weight), length = 10000)
  par(mfrow=c(1,2))
  qqplot(FXinv(u = ppoints(siteinput$n), b = PLB.bMLE.site.b[i], xmin = min(sitex),
                             xmax = max(sitex)), x , xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(A)", sites[i], ": Bounded power law Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
  qqplot(qlnorm(p = ppoints(siteinput$n), meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog), x, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(B)", sites[i], ": Log-normal Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))

}
dev.off()
