rm(list = ls())
library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(tidyverse)

source("coralsizedistfuncs.R")

#load data
stream <- read.csv("Pomeranz_estimated_dw.csv") #dw are all in grams

sites <- unique(stream$site)
nsites <- length(sites)

# create empty lists to store things in 
PLB.bMLE.site.b <- numeric(nsites)
PLB.minLL.site.b <- numeric(nsites)
siteb_plot <- list()
AICdf <- data.frame(site = sites, Number = NA, llBPL = NA, lllognorm = NA, AICBPL = NA, AIClognorm = NA)
gof <- data.frame(site = sites, lognormX2 = NA, lognormdf = NA, lognormP = NA, BPLX2 = NA, BPLdf = NA, BPLP = NA)
confintervals <- data.frame(site = sites, site_b = NA, lwr_CI95 = NA, upr_CI95 = NA, x_min = NA, x_max = NA, meanlog = NA, meanloglower95 = NA, meanlogupper95 = NA, sdlog = NA, sdloglower95 = NA, sdlogupper95 = NA)



for(i in 1:nsites){
  s <- sites[i]
  print(s)
  sitedata <- stream %>% filter(site == s)
  siteinput <- set.stream.params(sitedata$dw)
  print(summary(siteinput$biomass)) #WE CAN SEE THAT SAME DATA BEING USED EVERY TIME
  bML <- mle_b(Site == s, x = siteinput$biomass, sum_log_x = siteinput$sum.log.biomass, x_min = siteinput$min.biomass, x_max = siteinput$max.biomass)
  PLB.bMLE.site.b[i] <- bML[[1]]
  PLB.minLL.site.b[i] <- bML[[2]]$minimum
  thetalnorm <- estimatelognormal(x = sitedata$dw)

  #goodness-of-fit test for lognormal
  lngof <- lnormgof(x = sitedata$dw, mu = thetalnorm$meanlog, sigma = thetalnorm$sdlog)
  gof$lognormX2[i] <- round(lngof$X2,3)
  gof$lognormdf[i] <- lngof$df
  gof$lognormP[i] <- signif(lngof$P,3)
  
  #goodness-of-fit test for bounded power law
  bplgof <- BPLgof(x = sitedata$dw, b = bML[[1]], xmin = min(sitedata$dw), xmax = max(sitedata$dw))  #throwing error at "Portal"
  gof$BPLX2[i] <- round(bplgof$X2,3)
  gof$BPLdf[i] <- bplgof$df
  gof$BPLP[i] <- signif(bplgof$P,3)
  
  x <- sitedata$dw
  sitex = seq(min(sitedata$dw), max(sitedata$dw), length = 10000)
  
  #rank plot
  sitey.PLB <- (1 - pPLB(x = sitex, b = PLB.bMLE.site.b[i], xmin = min(sitex),
                         xmax = max(sitex))) * length(sitedata$dw)
  sitey.lnorm <- plnorm(q = sitex, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog, lower.tail = FALSE, log.p = FALSE) * length(sitedata$dw)

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
    geom_point(aes_(x = (sort(sitedata$dw, decreasing=TRUE)), y = (1:length(sitedata$dw))),
               color = "#666666", size = 2, alpha = 0.3) +
    scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
                       limits = c(0.25, max(table(stream$site)))) +
    scale_x_continuous(trans = 'log10',
                       limits = range(stream$dw))+
    geom_line(aes_(x = sitex, y = sitey.PLB), col = 'black', lwd = 1) +
    geom_line(aes_(x = sitex, y = sitey.lnorm), col = '#1B9E77', lwd = 1) +
    labs(tag = LETTERS[ ((i - 1) %% 4) + 1 ]) +
    annotate("text", x = 1e-05, y = 10, label = s) +
    annotate("text", x = 1e-05, y = 5, label = paste("italic(b)[PLB]==",(round(PLB.bMLE.site.b[i],2))), parse = T) +
       annotate("text", x = 1e-05, y = 2, label = paste("n =" ,(length(sitedata$dw)))) +
    theme_classic() +
    theme(axis.title = element_blank(), legend.position = "none")
  AICdf$Number[i] <- length(sitedata$dw)
  AICdf$lllognorm[i] <- lnormAIC(x)$lllognorm
  AICdf$AIClognorm[i] <- lnormAIC(x)$AIClognorm
  AICdf$llBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.biomass, xmax = siteinput$max.biomass, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$llBPL
  AICdf$AICBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.biomass, xmax = siteinput$max.biomass, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$AICBPL
  confintervals$site_b[i] <- round(PLB.bMLE.site.b[i],3)
  confintervals$lwr_CI95[i] <- round(slope.conf.int(PLB.bMLE.site.b[[i]], PLB.minLL.site.b[[i]], siteinput)[1],3)
  confintervals$upr_CI95[i] <- round(slope.conf.int(PLB.bMLE.site.b[[i]], PLB.minLL.site.b[[i]], siteinput)[2],3)
  confintervals$x_min[i] <- siteinput$min.biomass
  confintervals$x_max[i] <- round(siteinput$max.biomass, 3)
  confintervals$meanlog[i] <- round(thetalnorm$meanlog, 3)
  confintervals$meanloglower95[i] <- round(thetalnorm$meanlower, 3)
  confintervals$meanlogupper95[i] <- round(thetalnorm$meanupper, 3)
  confintervals$sdlog[i] <- round(thetalnorm$sdlog, 3)
  confintervals$sdloglower95[i] <- round(thetalnorm$sdlower, 3)
  confintervals$sdlogupper95[i] <- round(thetalnorm$sdupper, 3)
  
  }  


leftlabel <- grid::textGrob(expression(paste("Number of stream macroinverts with sizes", " ">=" ", italic("x"), "    ")), rot = 90)
bottomlabel <- grid::textGrob(expression(paste("Stream macroinverts biomass, ", italic("x"), ~(g))))

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
  "streambiomass.pdf",
  finalplot,
  width = 12,
  height = 8.5,
  units = "in"
)

write.csv(AICdf,'stream_AIC.csv')
write.csv(gof, 'gofstream.csv')
write.csv(confintervals, 'streamparam.csv')


# qqplots saved as pdf
pdf("qqplots_stream.pdf", width = 12, height = 6)
for(i in 1:nsites){
  s <- sites[i]
  sitedata <- stream %>% filter(site == s)
  siteinput <- set.stream.params(sitedata$dw)
  x <- sitedata$dw
  sitex = seq(min(sitedata$dw), max(sitedata$dw), length = 10000)
  par(mfrow=c(1,2))
  qqplot(FXinv(u = ppoints(siteinput$n), b = PLB.bMLE.site.b[i], xmin = min(sitex),
                             xmax = max(sitex)), x , xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(A)", sites[i], ": Bounded power law"), pch = 16, col = adjustcolor("black", 0.25))
  qqplot(qlnorm(p = ppoints(siteinput$n), meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog), x, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(B)", sites[i], ": Log-normal"), pch = 16, col = adjustcolor("black", 0.25))
 
}
dev.off()

