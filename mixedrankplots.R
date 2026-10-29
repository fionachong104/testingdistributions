rm(list = ls())
library(dplyr)
library(ggplot2)
library(gridExtra)
library(cowplot)

source("coralsizedistfuncs.R")

#load data
oneyeardf <- read.csv("oneyeardf.csv",row.names=1)
oneyeardf <- oneyeardf[is.na(oneyeardf$ROI.LabelCode), ] # to remove corals that were out of frame 

sites <- unique(oneyeardf$Site)[15]
nsites <- length(sites)

PLB.bMLE.site.b <- numeric(nsites)
PLB.minLL.site.b <- numeric(nsites)

siteb_plot_coral <- list()

for(i in 1:nsites){
  
  s <- sites[i]
  print(s)
  sitedata <- oneyeardf %>% filter(Site == s)
  siteinput <- set.params(sitedata$Area)
  bML <- mle_b(Site == s, x = siteinput$Area, sum_log_x = siteinput$sum.log.Area, x_min = siteinput$min.Area, x_max = siteinput$max.Area)
  PLB.bMLE.site.b[i] <- bML[[1]]
  PLB.minLL.site.b[i] <- bML[[2]]$minimum
  thetalnorm <- estimatelognormal(x = sitedata$Area)

  x <- sitedata$Area
  sitex = seq(min(sitedata$Area), max(sitedata$Area), length = 10000)
  
  #rank plot
  sitey.PLB <- (1 - pPLB(x = sitex, b = PLB.bMLE.site.b[i], xmin = min(sitex),
                         xmax = max(sitex))) * length(sitedata$Area)
  sitey.lnorm <- plnorm(q = sitex, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog, lower.tail = FALSE, log.p = FALSE) * length(sitedata$Area)

  model_df <- data.frame(
    sitex = sitex,
    value = c(sitey.PLB, sitey.lnorm),
    model = rep(c("Bounded Power Law", "Log-normal"), each = length(sitex))) 
  
  siteb_plot_coral[[i]] <- ggplot() +
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
    labs(tag = "A", y = expression("Number of corals with areas " >= italic(x)), x = expression(paste("Coral area, ", italic("x"), ~(cm^2))))+
    annotate("text", x = 10, y = 10, label = s) +
    annotate("text", x = 10, y = 5, label = paste("n =" ,(length(sitedata$Area)))) +
     theme_classic() +
     theme(legend.position = "none")
}  

tree <- read.csv("Ausplots_large_tree_survey_data.csv")
tree <- tree %>% drop_na(aboveGroundLiveBiomass_kilograms)  # to remove dead trees

sites <- unique(tree$siteName)[8]
nsites <- length(sites)

PLB.bMLE.site.b_tree <- numeric(nsites)
siteb_plot_tree <- list()

for(i in 1:nsites){
  s <- sites[i]
  sitedata <- tree %>% filter(siteName == s)
  siteinput <- set.fish.params(sitedata$aboveGroundLiveBiomass_kilograms)
  bML <- mle_b(Site == s, x = siteinput$biomass, sum_log_x = siteinput$sum.log.biomass,
               x_min = siteinput$min.biomass, x_max = siteinput$max.biomass)
  PLB.bMLE.site.b_tree[i] <- bML[[1]] 
  thetalnorm <- estimatelognormal(x = sitedata$aboveGroundLiveBiomass_kilograms)
  
  x <- sitedata$aboveGroundLiveBiomass_kilograms
  sitex = seq(min(sitedata$aboveGroundLiveBiomass_kilograms), max(sitedata$aboveGroundLiveBiomass_kilograms), length = 10000)
  par(mfrow=c(1,2))
  #rank plot   
  sitey.PLB <- (1 - pPLB(x = sitex, b = PLB.bMLE.site.b_tree[i], xmin = min(sitex),
                         xmax = max(sitex))) * length(sitedata$aboveGroundLiveBiomass_kilograms)
  sitey.lnorm <- plnorm(q = sitex, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog, lower.tail = FALSE, log.p = FALSE) * length(sitedata$aboveGroundLiveBiomass_kilograms)
  
  siteb_plot_tree[[i]] <- ggplot() +
    geom_point(aes_(x = (sort(sitedata$aboveGroundLiveBiomass_kilograms, decreasing=TRUE)), y = (1:length(sitedata$aboveGroundLiveBiomass_kilograms))),
               color = "#666666", size = 2, alpha = 0.3) +
    scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
                       limits = c(0.25, max(table(tree$siteName)))) +
    scale_x_continuous(trans = 'log10',
                       limits = range(tree$aboveGroundLiveBiomass_kilograms))+
    geom_line(aes_(x = sitex, y = sitey.PLB), col = 'black', lwd = 1) +
    geom_line(aes_(x = sitex, y = sitey.lnorm), col = '#1B9E77', lwd = 1) +
    labs(tag = "B", y = expression("Number of trees with biomass " >= italic(x))
, x = expression(paste("Tree biomass, ", italic("x"), ~(kg)))) +
    annotate("text", x = 100, y = 10, label = s) +
    annotate("text", x = 100, y = 5, label = paste("n =" ,(length(sitedata$aboveGroundLiveBiomass_kilograms)))) +
    theme_classic()
}

stream <- read.csv("Pomeranz_estimated_dw.csv") #dw are all in grams

sites <- unique(stream$site)[12]
nsites <- length(sites)

# create empty lists to store things in 
PLB.bMLE.site.b_stream <- numeric(nsites)
siteb_plot_stream <- list()

for(i in 1:nsites){
  s <- sites[i]
  print(s)
  sitedata <- stream %>% filter(site == s)
  siteinput <- set.stream.params(stream$dw)
  print(summary(siteinput$biomass)) #WE CAN SEE THAT SAME DATA BEING USED EVERY TIME
  bML <- mle_b(Site == s, x = siteinput$biomass, sum_log_x = siteinput$sum.log.biomass, x_min = siteinput$min.biomass, x_max = siteinput$max.biomass)
  PLB.bMLE.site.b_stream[i] <- bML[[1]] 
  thetalnorm <- estimatelognormal(x = sitedata$dw)

  x <- sitedata$dw
  sitex = seq(min(sitedata$dw), max(sitedata$dw), length = 10000)
  
  #rank plot
  sitey.PLB <- (1 - pPLB(x = sitex, b = PLB.bMLE.site.b[i], xmin = min(sitex),
                         xmax = max(sitex))) * length(sitedata$dw)
  sitey.lnorm <- plnorm(q = sitex, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog, lower.tail = FALSE, log.p = FALSE) * length(sitedata$dw)
  siteb_plot_stream[[i]] <- ggplot() +
    geom_point(aes_(x = (sort(sitedata$dw, decreasing=TRUE)), y = (1:length(sitedata$dw))),
               color = "#666666", size = 2, alpha = 0.3) +
    scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
                       limits = c(0.25, max(table(stream$site)))) +
    scale_x_continuous(trans = 'log10',
                       limits = range(stream$dw))+
    geom_line(aes_(x = sitex, y = sitey.PLB), col = 'black', lwd = 1) +
    geom_line(aes_(x = sitex, y = sitey.lnorm), col = '#1B9E77', lwd = 1) +
    labs(tag = "C", y = expression("Number of macroinverts with biomass " >= italic(x)), x = expression(paste("Stream macroinverts biomass, ", italic("x"), ~(g)))) +
    annotate("text", x = 1e-05, y = 10, label = s) +
    annotate("text", x = 1e-05, y = 5, label = paste("n =" ,(length(sitedata$dw)))) +
    theme_classic() 
  }
# +
#     theme(axis.title = element_blank())
# }  






# getting the legend from the first siteb_plot
legend_plot <- get_legend(
  siteb_plot_coral[[1]] +
    theme(legend.position = "top",
          legend.direction = "horizontal",
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16))
)

# putting the plots into a 5x4 grid

panel_grid <- arrangeGrob(
  grobs = list(siteb_plot_coral[[i]],siteb_plot_tree[[i]], siteb_plot_stream[[i]]),
  nrow  = 3,
  ncol  = 1,
)

final_plot <- plot_grid(
  legend_plot,
  panel_grid,
  ncol = 1,
  rel_heights = c(0.05, 1)
)

ggsave(
  filename = "mixedranks.pdf",
  final_plot,
  width = 8, height = 10, units = "in"
)

ggsave(
  filename = "mixedranks.png",
  final_plot,
  width = 8, height = 10, units = "in"
)

