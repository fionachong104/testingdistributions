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

AICdf <- data.frame(site = sites, llBPL = NA, llMSBPL = NA, lllognorm = NA, llmslognorm = NA, AICBPL = NA, AICMSBPL = NA, AIClognorm = NA, AICmslognorm = NA)

w <- 3648/35
v <- 2736/35

for(i in 1:nsites){
  s <- sites[i]
  sitedata <- oneyeardf %>% filter(Site == s)
  siteinput <- set.params(sitedata$Area)
  bML <- mle_b(Site == s, x = siteinput$Area, log_x = sitedata$log.Area, sum_log_x = siteinput$sum.log.Area,
                         x_min = siteinput$min.Area, x_max = siteinput$max.Area)
  PLB.bMLE.site.b[i] <- bML[[1]] 
  msbplfit <- estimatebMSBPL(x = sitedata$Area, w = w, v = v)
  MSPLB.bMLE.site.b[i] <- msbplfit$minimum
  thetalnorm <- estimatelognormal(x = sitedata$Area)
  thetaMSlnorm <- estimateMSlnorm(x = sitedata$Area, w = w, v = v)
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
  qqplot(qlnorm(p = ppoints(siteinput$n), meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog), x, xlab = "theoretical quantiles", ylab = "sample quantiles", main = paste("(C)", sites[i], ": Log-normal Q-Q plot"))
  qqline(x, distribution = function(p){qlnorm(p, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog)})
  qqplot(FMSlnorminv(u = ppoints(siteinput$n), mu = thetaMSlnorm$par[1], sigma = thetaMSlnorm$par[2], w = w, v = v), sitedata$Area, xlab = "theoretical quantiles", ylab = "sample quantiles", main = paste("(D)", sites[i], ": Minus-sampled log-normal Q-Q plot"))
  qqline(sitedata$Area, distribution = function(p){
  FMSlnorminv(p, mu = thetaMSlnorm$par[1], sigma = thetaMSlnorm$par[2], w = w, v = v)
  })

#rank plot   
  sitey.PLB <- (1 - pPLB(x = sitex, b = PLB.bMLE.site.b[i], xmin = min(sitex),
                       xmax = max(sitex))) * length(sitedata$Area)
  sitey.MSBPL <-  (1 - FMSBPL(x = sitex, b = MSPLB.bMLE.site.b[i], xmin = min(sitex),
                            w = w, v = v)) * length(sitedata$Area)
  sitey.MSlnorm <- (1-sapply(sitex, FUN = FMSlnorm, mu = thetaMSlnorm$par[1], sigma = thetaMSlnorm$par[2], w = w, v = v))* length(sitedata$Area)
  sitey.lnorm <- plnorm(q = sitex, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog, lower.tail = FALSE, log.p = FALSE) * length(sitedata$Area)
 # 1 - sapply(sitex.MSlnorm, FUN = FMSlnorm, mu = mu, sigma = sigma, w = w, v = v
  siteb_plot[[i]] <- ggplot() +
    geom_point(aes_(x = (sort(sitedata$Area, decreasing=TRUE)), y = (1:length(sitedata$Area))),
               color = "cadetblue", size = 2, alpha = 0.3) +
    scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000), 
                       limits = c(0.25, max(table(oneyeardf$Site)))) +
    scale_x_continuous(trans = 'log10', breaks = c(0,1,5,10,100,1000,10000),
                       limits = range(oneyeardf$Area))+
    geom_line(aes_(x = sitex, y = sitey.PLB), col = 'black', lwd = 1) +
    geom_line(aes_(x = sitex, y = sitey.MSBPL), col = 'red', lwd = 1) +
    geom_line(aes_(x = sitex, y = sitey.MSlnorm), col = 'blue', lwd = 1) +
    geom_line(aes_(x = sitex, y = sitey.lnorm), col = 'green', lwd = 1) +
    labs(tag = LETTERS[i]) +
    annotate("text", x = 10, y = 10, label = s) +
    annotate("text", x = 10, y = 3, label = bquote(paste(italic(b)[PLB]==.(round(PLB.bMLE.site.b[i],2))))) +
    annotate("text", x = 10, y = 1, label = bquote(paste(italic(b)[MSBPL]==.(round(MSPLB.bMLE.site.b[i],2))))) +
    annotate("text", x = 650, y = 800, label = bquote(n == .(length(sitedata$Area)))) +
    theme_classic() + 
    theme(axis.title = element_blank())
  AICdf$lllognorm[i] <- lnormAIC(x)$lllognorm
  AICdf$AIClognorm[i] <- lnormAIC(x)$AIClognorm
  AICdf$llBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.Area, xmax = siteinput$max.Area, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$llBPL
  AICdf$AICBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.Area, xmax = siteinput$max.Area, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$AICBPL
  AICdf$llmslognorm[i] <- MSlnormAIC(thetaMSlnorm)$llmslognorm
  AICdf$AICmslognorm[i] <- MSlnormAIC(thetaMSlnorm)$AICmslognorm
  AICdf$llMSBPL[i] <- MSBPLAIC(msbplfit)$llMSBPL
  AICdf$AICMSBPL[i] <- MSBPLAIC(msbplfit)$AICMSBPL
}

leftlabel <- grid::textGrob(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    ")), rot = 90)
bottomlabel <- grid::textGrob(expression(paste("Colony area, ", italic("x"), ~(cm^2))))

siteb_plot <- grid.arrange(grobs = siteb_plot, ncol = 4, 
             left = leftlabel,
             bottom = bottomlabel)
#ggsave(file = "siteb_plot.svg", plot = siteb_plot, width = 13, height = 9)

# site-wise SFD of log coral area
par(
  mfrow = c(5,4),
  mar = c(1,1.5,1,0),
  oma = c(4,4,2,2)
)  
for(i in 1:nsites){
  s <- sites[i]
  sitedata <- oneyeardf %>% filter(Site == s)
  custom_breaks <- c(-1:10)
  siteproportions <- hist(log(sitedata$Area), plot = FALSE, breaks = custom_breaks)
  plot(siteproportions, freq = FALSE, col = "darkgrey", main = "", xlab = "", ylab = "", axes = FALSE,  ylim = c(0,0.4), xlim = c(0,10))
  title(paste("(",LETTERS[i],")"," ",sites[i], sep= ""), cex.main = 1.5, line = -1.5)
  lines(density(log(sitedata$Area)), na.rm=TRUE, col = "blue", lty = "dashed", lwd = 2)
  abline(v=mean(log(sitedata$Area)), col = "red", lwd = 2, lty = "dashed")
  legend("topright", inset = .05, bty = "n", cex = 1.5, legend = bquote(n == .(length(sitedata$Area))))
  axis(side = 1, at=c(0,5,10), labels = c(0,5,10), cex.axis = 1.5)
  axis(side = 2, at=c(0,0.2,0.4), labels = c(0,0.2,0.4), cex.axis = 1.5)
  #xlab = expression(paste("Log coral area"~(cm^2))))
}
mtext(expression(paste("log(coral area/"*cm^2*")")) , side=1,line=3,outer=TRUE,cex=1.3)
mtext("Proportion of coral colonies", side=2,line=2,outer=TRUE,cex=1.3,las=0)

# 
# # saves the temp images from the plotting envrionment
# plots.dir.path <- list.files(tempdir(), pattern = "rs-graphics", full.names = TRUE)
# plots.png.paths <- list.files(plots.dir.path, pattern = ".png", full.names = TRUE)
# file.copy(from = plots.png.paths, to = "C:/Users/624225/OneDrive - hull.ac.uk/_BoxData/PhD/testingdistributions")


#write.csv(AICdf,'AIC.csv') 
