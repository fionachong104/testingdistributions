rm(list = ls())
library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)

source("coralsizedistfuncs.R")

#load data
oneyeardf <- read.csv("oneyeardf.csv",row.names=1)
oneyeardf <- oneyeardf[is.na(oneyeardf$ROI.LabelCode), ] # to remove corals that were out of frame 
axisscores <- read.csv("axisscores.csv", row.names=1)
axisscores <- axisscores[order(axisscores$PC1), ]

sites <- unique(oneyeardf$Site)
sites <- sites[match(row.names(axisscores), sites)] #re-ordering based on increasing PC1 scores
nsites <- length(sites)

PLB.bMLE.site.b <- numeric(nsites)
PLB.minLL.site.b <- numeric(nsites)

lnormparams <- data.frame(meanlog = numeric(nsites), sdlog = numeric(nsites))

siteb_plot <- list()
bplqq <- list()

AICdf <- data.frame(site = sites, Number = NA, llBPL = NA, lllognorm = NA, AICBPL = NA,  AIClognorm = NA)

gof <- data.frame(site = sites, lognormX2 = NA, lognormdf = NA, lognormP = NA, BPLX2 = NA, BPLdf = NA, BPLP = NA)


confintervals <- data.frame(site = sites, site_b = NA, lwr_CI95 = NA, upr_CI95 = NA, x_min = NA, x_max = NA)
confintervalslognormal <- data.frame(site = sites, meanlog = NA, meanloglower95 = NA, meanlogupper95 = NA, sdlog = NA, sdloglower95 = NA, sdlogupper95 = NA)
# 
# w <- 3648/35
# v <- 2736/35

for(i in 1:nsites){
  s <- sites[i]
  print(s)
  sitedata <- oneyeardf %>% filter(Site == s)
  siteinput <- set.params(sitedata$Area)
  bML <- mle_b(Site == s, x = siteinput$Area, sum_log_x = siteinput$sum.log.Area, x_min = siteinput$min.Area, x_max = siteinput$max.Area)
  PLB.bMLE.site.b[i] <- bML[[1]]
  PLB.minLL.site.b[i] <- bML[[2]]$minimum

   thetalnorm <- estimatelognormal(x = sitedata$Area)
   lnormparams$meanlog[i] <- thetalnorm$meanlog
   lnormparams$sdlog[i] <- thetalnorm$sdlog
   
  #goodness-of-fit test for lognormal
  lngof <- lnormgof(x = sitedata$Area, mu = thetalnorm$meanlog, sigma = thetalnorm$sdlog)
  gof$lognormX2[i] <- round(lngof$X2,3)
  gof$lognormdf[i] <- lngof$df
  gof$lognormP[i] <- signif(lngof$P,3)
  
  #goodness-of-fit test for bounded power law
  bplgof <- BPLgof(x = sitedata$Area, b = bML[[1]], xmin = min(sitedata$Area), xmax = max(sitedata$Area))
  gof$BPLX2[i] <- signif(bplgof$X2,3)
  gof$BPLdf[i] <- bplgof$df
  gof$BPLP[i] <- signif(bplgof$P,3)
  
    
  x <- sitedata$Area
  sitex = seq(min(sitedata$Area), max(sitedata$Area), length = 10000)

#rank plot
  sitey.PLB <- (1 - pPLB(x = sitex, b = PLB.bMLE.site.b[i], xmin = min(sitex),
                       xmax = max(sitex))) * length(sitedata$Area)
  sitey.lnorm <- plnorm(q = sitex, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog, lower.tail = FALSE, log.p = FALSE) * length(sitedata$Area)
   model_df <- data.frame(
    sitex = sitex,
    value = c(sitey.PLB, sitey.lnorm),
    model = rep(c("Bounded Power Law", "Log-normal"), each = length(sitex)) 
                
  )

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
    geom_point(aes_(x = (sort(sitedata$Area, decreasing=TRUE)), y = (1:length(sitedata$Area))),
               color = "#666666", size = 2, alpha = 0.2) +
    scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
                       limits = c(0.25, max(table(oneyeardf$Site)))) +
    scale_x_continuous(trans = 'log10', breaks = c(0,1,5,10,100,1000,10000),
                       limits = range(oneyeardf$Area))+
    labs(tag = LETTERS[i])+ 
    annotate("text", x = 10, y = 10, label = s) +
    annotate("text", x = 10, y = 5, label = paste("n =" ,(length(sitedata$Area)))) +
    theme_classic() +
    theme(axis.title = element_blank(), legend.position = "none")
  AICdf$Number[i] <- length(sitedata$Area)
  AICdf$lllognorm[i] <- round(lnormAIC(x)$lllognorm,3)
  AICdf$AIClognorm[i] <- round(lnormAIC(x)$AIClognorm,3)
  AICdf$llBPL[i] <- round(BPLAIC(C = getC(xmin = siteinput$min.Area, xmax = siteinput$max.Area, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$llBPL,3)
  AICdf$AICBPL[i] <- round(BPLAIC(C = getC(xmin = siteinput$min.Area, xmax = siteinput$max.Area, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$AICBPL,3)
  confintervals$site_b[i] <- round(PLB.bMLE.site.b[i],3)
  confintervals$lwr_CI95[i] <- round(slope.conf.int.coral(PLB.bMLE.site.b[[i]], PLB.minLL.site.b[[i]], siteinput)[1],3)
  confintervals$upr_CI95[i] <- round(slope.conf.int.coral(PLB.bMLE.site.b[[i]], PLB.minLL.site.b[[i]], siteinput)[2],3)
  confintervals$x_min[i] <- round(siteinput$min.Area, 3)
  confintervals$x_max[i] <- round(siteinput$max.Area, 3)
  confintervalslognormal$meanlog[i] <- round(thetalnorm$meanlog, 3)
  confintervalslognormal$meanloglower95[i] <- round(thetalnorm$meanlower, 3)
  confintervalslognormal$meanlogupper95[i] <- round(thetalnorm$meanupper, 3)
  confintervalslognormal$sdlog[i] <- round(thetalnorm$sdlog, 3)
  confintervalslognormal$sdloglower95[i] <- round(thetalnorm$sdlower, 3)
  confintervalslognormal$sdlogupper95[i] <- round(thetalnorm$sdupper, 3)
}  


# getting the legend from the first siteb_plot
legend_plot <- get_legend(
  siteb_plot[[1]] +
    theme(legend.position = "top",
          legend.direction = "horizontal",
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16))
)

# putting the plots into a 5x4 grid

leftlabel <- grid::textGrob(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    ")), rot = 90,
                            gp = grid::gpar(fontsize = 16))
bottomlabel <- grid::textGrob(expression(paste("Coral area, ", italic("x"), ~(cm^2))), gp = grid::gpar(fontsize = 16))

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
  filename = "coralarea.pdf", 
  final_plot, 
  width = 15, height = 13
)


parameters <- merge(confintervals, confintervalslognormal)

write.csv(AICdf,'coralAIC.csv') 
write.csv(gof,'coralgof.csv') 
write.csv(parameters, 'coralparam.csv')


# qqplots saved as pdf
pdf("qqplots_corals.pdf", width = 12, height = 6)
for(i in 1:nsites){
  s <- sites[i]
  print(s)
  sitedata <- oneyeardf %>% filter(Site == s)
  siteinput <- set.params(sitedata$Area)
  x <- sitedata$Area
  sitex = seq(min(sitedata$Area), max(sitedata$Area), length = 10000)
  par(mfrow=c(1,2))
  bplqq[[i]] <- qqplot(FXinv(u = ppoints(siteinput$n), b = PLB.bMLE.site.b[i], xmin = min(sitex),
                             xmax = max(sitex)), x , xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(A)", sites[i], ": Bounded power law"), pch = 16, col = adjustcolor("black", 0.25))
 qqplot(qlnorm(p = ppoints(siteinput$n), meanlog = lnormparams$meanlog[i], sdlog = lnormparams$sdlog[i]), x, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(B)", sites[i], ": Log-normal"), pch = 16, col = adjustcolor("black", 0.25))

}
dev.off()


