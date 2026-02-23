rm(list = ls())
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)

source("coralsizedistfuncs.R")

#load data
stream <- read.csv("Pomeranz_estimated_dw.csv")
sites <- unique(stream$site)
#sites <- sites[match(row.names(axisscores), sites)] #re-ordering based on increasing PC1 scores
nsites <- length(sites)


# create empty lists to store things in 
PLB.bMLE.site.b <- numeric(nsites)
siteb_plot <- list()
bplqq <- list()
AICdf <- data.frame(site = sites, Number = NA, llBPL = NA, lllognorm = NA, AICBPL = NA, AIClognorm = NA)


for(i in 1:nsites){
  s <- sites[i]
  sitedata <- stream %>% filter(site == s)
  siteinput <- set.stream.params(sitedata$dw)
  bML <- mle_b(Site == s, x = siteinput$biomass, sum_log_x = siteinput$sum.log.biomass,
               x_min = siteinput$min.biomass, x_max = siteinput$max.biomass)
  PLB.bMLE.site.b[i] <- bML[[1]] 
  thetalnorm <- estimatelognormal(x = sitedata$dw)
  x <- sitedata$dw
  sitex = seq(min(sitedata$dw), max(sitedata$dw), length = 10000)
  par(mfrow=c(1,2))
  bplqq[[i]] <- qqplot(FXinv(u = ppoints(siteinput$n), b = PLB.bMLE.site.b[i], xmin = min(sitex),
                             xmax = max(sitex)), x , xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(A)", sites[i], ": Power law Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
  #qqline(x, distribution = function(p){
  #  FXinv(p, b = PLB.bMLE.site.b[i], xmin = min(sitex), xmax = max(sitex))
  #}, untf=T)
  qqplot(qlnorm(p = ppoints(siteinput$n), meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog), x, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(B)", sites[i], ": Log-normal Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
  #  qqline(x, distribution = function(p){qlnorm(p, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog)}, untf=T)
  #rank plot   
  sitey.PLB <- (1 - pPLB(x = sitex, b = PLB.bMLE.site.b[i], xmin = min(sitex),
                         xmax = max(sitex))) * length(sitedata$dw)
  sitey.lnorm <- plnorm(q = sitex, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog, lower.tail = FALSE, log.p = FALSE) * length(sitedata$dw)
  siteb_plot[[i]] <- ggplot() +
    geom_point(aes_(x = (sort(sitedata$dw, decreasing=TRUE)), y = (1:length(sitedata$dw))),
               color = "#666666", size = 2, alpha = 0.3) +
    scale_y_continuous(trans = 'log10', breaks = c(1,10,100,300, 500), 
                       limits = c(0.25, max(table(stream$site)))) +
    scale_x_continuous(trans = 'log10',#breaks = c(-0.00000001, -0.000001, -0.001, -0.1),
                       limits = range(stream$dw))+
    geom_line(aes_(x = sitex, y = sitey.PLB), col = 'black', lwd = 1) +
    geom_line(aes_(x = sitex, y = sitey.lnorm), col = '#1B9E77', lwd = 1) +
    labs(tag = LETTERS[i]) +
    annotate("text", x =0.00001, y = 1, label = s) +
    annotate("text", x = 0.00001, y = 10, label = paste("italic(b)[PLB]==",(round(PLB.bMLE.site.b[i],2))), parse = T) +
    annotate("text", x = 0.00001, y =50, label = paste("n =" ,(length(sitedata$dw)))) +
    theme_classic() + 
    theme(axis.title = element_blank())
  AICdf$Number[i] <- length(sitedata$dw)
  AICdf$lllognorm[i] <- lnormAIC(x)$lllognorm
  AICdf$AIClognorm[i] <- lnormAIC(x)$AIClognorm
  AICdf$llBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.biomass, xmax = siteinput$max.biomass, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$llBPL
  AICdf$AICBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.biomass, xmax = siteinput$max.biomass, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$AICBPL
}

leftlabel <- grid::textGrob(expression(paste("Number of macroinvertebrates", " ">=" ", italic("x"), "    ")), rot = 90)
bottomlabel <- grid::textGrob(expression(paste("Stream macroinvertebrate dry weight, ", italic("x"), ~(g))))


ggsave(
  filename = "streaminvertbiomass.pdf", 
  plot = marrangeGrob(grobs= siteb_plot, nrow=2, ncol=2,
                      left = leftlabel,
                      bottom = bottomlabel,
                      layout_matrix = rbind(c(1,2), c(3,4))), 
  width = 15, height = 9
)


#repeat using lengths

# create empty lists to store things in 
PLB.bMLE.site.b_length <- numeric(nsites)
siteb_plot_length <- list()
bplqq_length <- list()
AICdf_length <- data.frame(site = sites, Number = NA, llBPL = NA, lllognorm = NA, AICBPL = NA, AIClognorm = NA)


for(i in 1:nsites){
  s <- sites[i]
  sitedata <- stream %>% filter(site == s)
  siteinput_length <- set.stream.params(sitedata$linear_meas)
  bML <- mle_b(Site == s, x = siteinput_length$biomass, sum_log_x = siteinput_length$sum.log.biomass,
               x_min = siteinput_length$min.biomass, x_max = siteinput_length$max.biomass)
  PLB.bMLE.site.b_length[i] <- bML[[1]] 
  thetalnorm_length<- estimatelognormal(x = sitedata$linear_meas)
  x_length <- sitedata$linear_meas
  sitex_length = seq(min(sitedata$linear_meas), max(sitedata$linear_meas), length = 10000)
  par(mfrow=c(1,2))
  bplqq[[i]] <- qqplot(FXinv(u = ppoints(siteinput_length$n), b = PLB.bMLE.site.b[i], xmin = min(sitex_length),
                             xmax = max(sitex_length)), x_length , xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(A)", sites[i], ": Power law Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
  #qqline(x, distribution = function(p){
  #  FXinv(p, b = PLB.bMLE.site.b[i], xmin = min(sitex), xmax = max(sitex))
  #}, untf=T)
  qqplot(qlnorm(p = ppoints(siteinput_length$n), meanlog = thetalnorm_length$meanlog, sdlog = thetalnorm_length$sdlog), x_length, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(B)", sites[i], ": Log-normal Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
  #  qqline(x, distribution = function(p){qlnorm(p, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog)}, untf=T)
  #rank plot   
  sitey.PLB_length <- (1 - pPLB(x = sitex_length, b = PLB.bMLE.site.b[i], xmin = min(sitex_length),
                         xmax = max(sitex_length))) * length(sitedata$linear_meas)
  sitey.lnorm_length <- plnorm(q = sitex_length, meanlog = thetalnorm_length$meanlog, sdlog = thetalnorm_length$sdlog, lower.tail = FALSE, log.p = FALSE) * length(sitedata$linear_meas)
  siteb_plot_length[[i]] <- ggplot() +
    geom_point(aes_(x = (sort(sitedata$linear_meas, decreasing=TRUE)), y = (1:length(sitedata$linear_meas))),
               color = "#666666", size = 2, alpha = 0.3) +
    scale_y_continuous(trans = 'log10', breaks = c(1,10,100,300, 500), 
                       limits = c(0.25, max(table(stream$site)))) +
    scale_x_continuous(trans = 'log10',#breaks = c(-0.00000001, -0.000001, -0.001, -0.1),
                       limits = range(stream$linear_meas))+
    geom_line(aes_(x = sitex_length, y = sitey.PLB_length), col = 'black', lwd = 1) +
    geom_line(aes_(x = sitex_length, y = sitey.lnorm_length), col = '#1B9E77', lwd = 1) +
    labs(tag = LETTERS[i]) +
    annotate("text", x =1, y = 1, label = s) +
    annotate("text", x = 1, y = 10, label = paste("italic(b)[PLB]==",(round(PLB.bMLE.site.b_length[i],2))), parse = T) +
    annotate("text", x = 1, y =50, label = paste("n =" ,(length(sitedata$linear_meas)))) +
    theme_classic() + 
    theme(axis.title = element_blank())
  AICdf_length$Number[i] <- length(sitedata$linear_meas)
  AICdf_length$lllognorm[i] <- lnormAIC(x_length)$lllognorm
  AICdf_length$AIClognorm[i] <- lnormAIC(x_length)$AIClognorm
  AICdf_length$llBPL[i] <- BPLAIC(C = getC(xmin = siteinput_length$min.biomass, xmax = siteinput_length$max.biomass, b = PLB.bMLE.site.b_length[i]), b = PLB.bMLE.site.b_length[i], x = x_length)$llBPL
  AICdf_length$AICBPL[i] <- BPLAIC(C = getC(xmin = siteinput_length$min.biomass, xmax = siteinput_length$max.biomass, b = PLB.bMLE.site.b_length[i]), b = PLB.bMLE.site.b_length[i], x = x_length)$AICBPL
}

leftlabel_length <- grid::textGrob(expression(paste("Number of macroinvertebrates", " ">=" ", italic("x"), "    ")), rot = 90)
bottomlabel_length <- grid::textGrob(expression(paste("Stream macroinvertebrate lengths, ", italic("x"), ~(mm))))


ggsave(
  filename = "streaminvertlengths.pdf", 
  plot = marrangeGrob(grobs= siteb_plot_length, nrow=2, ncol=2,
                      left = leftlabel_length,
                      bottom = bottomlabel_length,
                      layout_matrix = rbind(c(1,2), c(3,4))), 
  width = 15, height = 9
)

write.csv(AICdf,'stream_AIC.csv')
write.csv(AICdf_length,'stream_length_AIC.csv')

# # # saves the temp images from the plotting envrionment
# plots.dir.path <- list.files(tempdir(), pattern = "rs-graphics", full.names = TRUE)
# plots.png.paths <- list.files(plots.dir.path, pattern = ".png", full.names = TRUE)
# file.copy(from = plots.png.paths, to = "C:/Users/624225/OneDrive - hull.ac.uk/_BoxData/PhD/testingdistributions/qqplots/fish_logged")
