rm(list = ls())
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(moments) 

source("coralsizedistfuncs.R")

library(readxl)
#fish_size_spectra_data <- read_excel("C:/Users/624225/OneDrive - hull.ac.uk/_BoxData/PhD/testingdistributions/Carvalho21_code/fish_size_spectra_data.xlsx")
fishbiomass <- read_delim("fishbiomass_oneyear.csv")
axisscores <- read.csv("axisscores.csv", row.names=1)
axisscores <- axisscores[order(axisscores$PC1), ]

#reformatting Maria's fish data
#starting with BR2010
# BR_2010 <- as_tibble(BR_2010)
# TN <- BR_2010 %>%
#   group_by(Fish) %>%
#   summarise(total = sum(Number)) # this calculates how many fish there are, by species
fishbiomass <- as_tibble(fishbiomass)
Check_fishbiomass <- fishbiomass %>%
  group_by(Fish) %>%
  summarise(total = sum(Number))
#create data frame with each individual fish taking a row.
##make note of whether it is inividual biomass or not - each biomass value should be of a single fish not with them all grouped
# BRdf <- BR_2010 %>% 
#   uncount(Number) %>% 
#   group_by(Fish, Size) 

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
    #CV = round((sd/mean_log_biomass)*100, digits = 3),  removing this as our emans a negative so quite misleading
    skewness = round(skewness(log(individual_biomass_kg)), digits = 3),
    kurtosis = round(kurtosis(log(individual_biomass_kg)), digits = 3))
all_fish<- all_fish[match(rownames(axisscores), all_fish$Site), ] # order by PC1 scores 
write.csv(all_fish,'fishsummarystats.csv') 

sites <- unique(allfish$Site) #carn vs herb
sites <- sites[match(row.names(axisscores), sites)] #re-ordering based on increasing PC1 scores
nsites <- length(sites)

# create empty lists to store things in 
PLB.bMLE.site.b <- numeric(nsites)
siteb_plot <- list()
bplqq <- list()
AICdf <- data.frame(site = sites, llBPL = NA, lllognorm = NA, AICBPL = NA, AIClognorm = NA)


for(i in 1:nsites){
  s <- sites[i]
  sitedata <- allfish %>% filter(Site == s)
  siteinput <- set.fish.params(sitedata$individual_biomass_kg)
  bML <- mle_b(Site == s, x = siteinput$biomass, sum_log_x = siteinput$sum.log.biomass,
               x_min = siteinput$min.biomass, x_max = siteinput$max.biomass)
  PLB.bMLE.site.b[i] <- bML[[1]] 
  thetalnorm <- estimatelognormal(x = sitedata$individual_biomass_kg)
  x <- sitedata$individual_biomass_kg
  sitex = seq(min(sitedata$individual_biomass_kg), max(sitedata$individual_biomass_kg), length = 10000)
  par(mfrow=c(1,2))
  bplqq[[i]] <- qqplot(FXinv(u = ppoints(siteinput$n), b = PLB.bMLE.site.b[i], xmin = min(sitex),
                             xmax = max(sitex)), x , xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(A)", sites[i], ": Power law Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
  qqline(x, distribution = function(p){
    FXinv(p, b = PLB.bMLE.site.b[i], xmin = min(sitex), xmax = max(sitex))
  }, untf=T)
  qqplot(qlnorm(p = ppoints(siteinput$n), meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog), x, xlab = "Theoretical Quantiles", ylab = "Sample Quantiles", log = "xy", main = paste("(B)", sites[i], ": Log-normal Q-Q plot"), pch = 16, col = adjustcolor("black", 0.25))
  qqline(x, distribution = function(p){qlnorm(p, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog)}, untf=T)
  #rank plot   
  sitey.PLB <- (1 - pPLB(x = sitex, b = PLB.bMLE.site.b[i], xmin = min(sitex),
                         xmax = max(sitex))) * length(sitedata$individual_biomass_kg)
  sitey.lnorm <- plnorm(q = sitex, meanlog = thetalnorm$meanlog, sdlog = thetalnorm$sdlog, lower.tail = FALSE, log.p = FALSE) * length(sitedata$individual_biomass_kg)
  siteb_plot[[i]] <- ggplot() +
    geom_point(aes_(x = (sort(sitedata$individual_biomass_kg, decreasing=TRUE)), y = (1:length(sitedata$individual_biomass_kg))),
               color = "#666666", size = 2, alpha = 0.3) +
    scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000), 
                       limits = c(0.25, max(table(allfish$Site)))) +
    scale_x_continuous(trans = 'log10',#breaks = c(-10,0,1,10,100),
                        limits = range(allfish$individual_biomass_kg))+
    geom_line(aes_(x = sitex, y = sitey.PLB), col = 'black', lwd = 1) +
    geom_line(aes_(x = sitex, y = sitey.lnorm), col = '#1B9E77', lwd = 1) +
    labs(tag = LETTERS[i]) +
    annotate("text", x = 0.002, y = 10, label = s) +
    annotate("text", x = 0.001, y = 3, label = bquote(paste(italic(b)[PLB]==.(round(PLB.bMLE.site.b[i],2))))) +
    annotate("text", x = 100, y =1000, label = bquote(n == .(length(sitedata$individual_biomass_kg)))) +
    theme_classic() + 
    theme(axis.title = element_blank())
  AICdf$lllognorm[i] <- lnormAIC(x)$lllognorm
  AICdf$AIClognorm[i] <- lnormAIC(x)$AIClognorm
  AICdf$llBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.biomass, xmax = siteinput$max.biomass, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$llBPL
  AICdf$AICBPL[i] <- BPLAIC(C = getC(xmin = siteinput$min.biomass, xmax = siteinput$max.biomass, b = PLB.bMLE.site.b[i]), b = PLB.bMLE.site.b[i], x = x)$AICBPL
}

leftlabel <- grid::textGrob(expression(paste("Number of fish with sizes", " ">=" ", italic("x"), "    ")), rot = 90)
#bottomlabel <- grid::textGrob(expression(paste("Fish biomass, ", italic("x"), ~(kg))))
bottomlabel <- grid::textGrob(expression(paste("Fish biomass, ", italic("x"), ~(kg))))

siteb_plot <- grid.arrange(grobs = siteb_plot, ncol = 4,
                           left = leftlabel,
                           bottom = bottomlabel)


write.csv(AICdf,'fish_AIC.csv')

par(
  mfrow = c(5,4),
  mar = c(1,1.5,1,0),
  oma = c(4,4,2,2)
)
for(i in 1:nsites){
  s <- sites[i]
  sitedata <- allfish %>% filter(Site == s)
  custom_breaks <- c(-14:3)
  siteproportions <- hist(log(sitedata$individual_biomass_kg), plot = FALSE) #breaks = custom_breaks)
  plot(siteproportions, freq = FALSE, col = "darkgrey", main = "", xlab = "", ylab = "", axes = FALSE,  ylim = c(0,0.9), xlim = c(-15,3))
  title(paste("(",LETTERS[i],")"," ",sites[i], sep= ""), cex.main = 1.5, line = -1.5)
  lines(density(log(sitedata$individual_biomass_kg)), na.rm=TRUE, col = "blue", lty = "dashed", lwd = 2)
  abline(v=mean(log(sitedata$individual_biomass_kg)), col = "red", lwd = 2, lty = "dashed")
  legend("topright", inset = .05, bty = "n", cex = 1.5, legend = bquote(n == .(length(sitedata$individual_biomass_kg))))
  axis(side = 1, at=c(-14,-12,-10,-8, -6, -4,-2, 0, 2), labels = c(14,-12,-10,-8, -6, -4,-2, 0, 2), cex.axis = 1.5)
  axis(side = 2, at=c(0,0.45,0.9), labels = c(0,0.45,0.9), cex.axis = 1.5)}
  #xlab = expression(paste("Log coral area"~(cm^2))))

mtext(expression(paste("log(fish biomass (kg))")) , side=1,line=3,outer=TRUE,cex=1.3)
mtext("Density", side=2,line=2,outer=TRUE,cex=1.3,las=0)


# # # saves the temp images from the plotting envrionment
# plots.dir.path <- list.files(tempdir(), pattern = "rs-graphics", full.names = TRUE)
# plots.png.paths <- list.files(plots.dir.path, pattern = ".png", full.names = TRUE)
# file.copy(from = plots.png.paths, to = "C:/Users/624225/OneDrive - hull.ac.uk/_BoxData/PhD/testingdistributions/qqplots/fish_logged")
