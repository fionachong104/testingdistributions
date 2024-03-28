rm(list = ls())
library(dplyr)
library(ggplot2)
library(gridExtra)

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
MSPLB.bMLE.site.b <- numeric(nsites)

siteb_plot <- list()
bplqq <- list()
msbplqq <-list()

lognormalAIC <- list()
boundedpowerlawAIC <- list()
mslognormalAIC <- list()

w <- 3648/35
v <- 2736/35

for(i in 1:nsites){
  s <- sites[i]
  sitedata <- oneyeardf %>% filter(Site == s)
  siteinput <- set.params(sitedata$Area)
  bML <- mle_b(Site == s, x = siteinput$Area, log_x = sitedata$log.Area, sum_log_x = siteinput$sum.log.Area,
                         x_min = siteinput$min.Area, x_max = siteinput$max.Area)
  PLB.bMLE.site.b[i] <- bML[[1]] 
  MSPLB.bMLE.site.b[i] <- estimatebMSBPL(x = sitedata$Area, w = w, v = v)$minimum
  thetaML <- estimateMSlnorm(x = sitedata$Area, w = w, v = v)
  x <- sitedata$Area
  sitex = seq(min(sitedata$Area), max(sitedata$Area), length = 1000)
  par(mfrow=c(2,2))
  bplqq[[i]] <- qqplot(FXinv(u = ppoints(siteinput$n), b = PLB.bMLE.site.b[i], xmin = min(sitex),
               xmax = max(sitex)), x , xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", main = paste("(A)", sites[i], ": Power law Q-Q plot"))
  qqline(x, distribution = function(p){
    FXinv(p, b = PLB.bMLE.site.b[i], xmin = min(sitex), xmax = max(sitex))
  })
  msbplqq[[i]] <- qqplot(FMSBPLinv(u = ppoints(siteinput$n), b = MSPLB.bMLE.site.b[i], xmin = siteinput$min.Area, w = w, v = v), x, xlab = "Theoretical quantiles", ylab = "Sample quantiles", main = paste("(B)", sites[i], ": Minus-sampled bounded power law Q-Q plot"))
  qqline(x, distribution = function(p){
    FMSBPLinv(p, b =  MSPLB.bMLE.site.b[i], xmin = siteinput$min.Area, w = w, v = v)
  })
  #hist(log(sitedata$Area), main = paste("(C)", sites[i], ": Size-frequency distribution"), xlab = expression(paste("Log coral area"~(cm^2))))
  qqplot(qlnorm(p = ppoints(siteinput$n), meanlog = mean(log(x)), sdlog = sd(log(x))), x, xlab = "theoretical quantiles", ylab = "sample quantiles", main = " lognormal Q-Q plot")
  qqline(x, distribution = function(p){qlnorm(p, meanlog = mean(log(x)), sdlog = sd(log(x)))})
  qqplot(FMSlnorminv(u = ppoints(siteinput$n), mu = thetaML$par[1], sigma = thetaML$par[2], w = w, v = v), sitedata$Area, xlab = "theoretical quantiles", ylab = "sample quantiles", main = "minus-sampled lognormal Q-Q plot")
  qqline(sitedata$Area, distribution = function(p){
    FMSlnorminv(p, mu = thetaML$par[1], sigma = thetaML$par[2], w = w, v = v)
  })
  sitey.PLB = (1 - pPLB(x = sitex, b = PLB.bMLE.site.b[i], xmin = min(sitex),
                       xmax = max(sitex))) * length(sitedata$Area)
  sitey.MSBPL = (1 - FMSBPL(x = sitex, b = MSPLB.bMLE.site.b[i], xmin = min(sitex),
                            w = w, v = v)) * length(sitedata$Area)
  siteb_plot[[i]] <- ggplot() +
    geom_point(aes_(x = (sort(sitedata$Area, decreasing=TRUE)), y = (1:length(sitedata$Area))),
               color = "cadetblue", size = 2, alpha = 0.3) +
    scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000), 
                       limits = c(0.25, max(table(oneyeardf$Site)))) +
    scale_x_continuous(trans = 'log10', breaks = c(0,1,5,10,100,1000,10000),
                       limits = range(oneyeardf$Area))+
    geom_line(aes_(x = sitex, y = sitey.PLB), col = 'black', lwd = 1) +
    geom_line(aes_(x = sitex, y = sitey.MSBPL), col = 'red', lwd = 1) +
    labs(tag = LETTERS[i]) +
    annotate("text", x = 10, y = 10, label = s) +
    annotate("text", x = 10, y = 3, label = bquote(paste(italic(b)[PLB]==.(round(PLB.bMLE.site.b[i],2))))) +
    annotate("text", x = 10, y = 1, label = bquote(paste(italic(b)[MSBPL]==.(round(MSPLB.bMLE.site.b[i],2))))) +
    theme_classic() + 
    theme(axis.title = element_blank())
  lognormalAIC[[i]] <- lnormAIC(x)
  boundedpowerlawAIC[[i]] <- BPLAIC(C = getC(xmin = siteinput$min.Area, xmax = siteinput$max.Area, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)
  mslognormalAIC[[i]] <- MSlnormAIC(thetaML = thetaML)
}

leftlabel <- grid::textGrob(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    ")), rot = 90)
bottomlabel <- grid::textGrob(expression(paste("Colony area, ", italic("x"), ~(cm^2))))

siteb_plot <- grid.arrange(grobs = siteb_plot, ncol = 4, 
             left = leftlabel,
             bottom = bottomlabel)

ggsave(file = "siteb_plot.svg", plot = siteb_plot, width = 13, height = 9)

# saves the temp images from the plotting envrionment
# plots.dir.path <- list.files(tempdir(), pattern = "rs-graphics", full.names = TRUE)
# plots.png.paths <- list.files(plots.dir.path, pattern = ".png", full.names = TRUE)
# file.copy(from = plots.png.paths, to = "C:/Users/624225/OneDrive - hull.ac.uk/_BoxData/PhD/testingdistributions")

# compare AIC easily 
boundedpowerlawAIC <- do.call(rbind.data.frame,boundedpowerlawAIC)
lognormalAIC <- do.call(rbind.data.frame,lognormalAIC)
mslognormalAIC <- do.call(rbind.data.frame, mslognormalAIC)
AICdf <- cbind(boundedpowerlawAIC, lognormalAIC, mslognormalAIC)
row.names(AICdf) <- sites
write.csv(AICdf,'AIC.csv') 
