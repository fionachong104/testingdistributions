rm(list = ls())
library(RColorBrewer)

source("coralsizedistfuncs.R")

#dimensions of sampling window
# w <- 5
# v <- 2 #ASSUME v <= w

# images dimensions from CANON camera (w x v = 3648 x 2736)
# 10 cm = roughly 350 pixels; 1 m = 3500; 1 cm = 35; so dividing by 35 to see what that translates to in cm x cm
w <- 3648/35
v <- 2736/35

#DATA FILES NEED TO BE IN THE PROJECT FOR THIS TO WORK
#load data
#oneyeardf <- read.csv("oneyeardf.csv")
#oneyeardf <- oneyeardf[is.na(oneyeardf$ROI.LabelCode), ] # subset by removing corals that are out of frame
#axisscores <- read.csv("axisscores.csv")
#axisscores <- axisscores[order(axisscores$PC1), ] # order by PC1 from low to high

#sample from the sites; loop to repeat for all sites
#subsample from the minus sampled distribution



#parameters of bounded power law
# need to make a function that estimates the parameters for each site
xmin <- 0.56 #1e-1 #min area FC: 0.56 cm^2 from oneyeardf
xmax <- 4569 #pi * 3^2 #max area FC: 4569 cm^2 from oneyeardf
xmaxminus <- pi / 4 * v^2 #max area of circle that can fit in window
b <- -1.7 #power law exponent
C <- getC(xmin = xmin, xmax = xmax, b = b) #normalization constant
n <- 1e4 #number of circles
nplot <- 100 #number to plot
x <- simulateboundedpowerlaw(n = n, b = b, xmin = xmin, xmax = xmax) #simulated areas

r <- radiusfromarea(x) #colony radii (assume circles)
mycolors <- brewer.pal(3, "Dark2")

alpha <- runif(n = n, min = 0, max = w) #(alpha, beta) uniform random points in sampling window
beta <- runif(n = n, min = 0, max = v)

par(mfrow = c(2, 3))
plot(c(0, w), c(0, v), type = "n", asp = 1, xlab = "", ylab = "")
plotframe(w = w, v = v)
isinframe <- inframe(alpha = alpha, beta = beta, r = r, w = w, v = v)
for(i in 1:nplot){#drawing all the circles can be slow so make nplot not too large
  plotcircle(alpha = alpha[i], beta = beta[i], r = r[i], isinframe = isinframe[i])
}

hist(x, freq = FALSE, main = "bounded power law", breaks = 100, xlab = expression(italic(x)), ylab = "density")
xseq <- seq(from = xmin, to = xmax, length.out = 1e4)
fx <- dboundedpowerlaw(x = xseq, b = b, C = C, xmin = xmin, xmax = xmax)
lines(xseq, fx, col = mycolors[1], lwd = 2, lty = "solid") #density curve for bounded power law
legend("topright", bty = "n", lwd = 2, lty = c("solid"), col = c(mycolors[1]), legend = c("calculated"))

hist(x[isinframe], freq = FALSE, main = "minus-sampled bounded power law", breaks = 100, xlab = expression(italic(x)), ylab = "density", ylim = c(0, 0.2)) #changing y-axis limits helps us see tail
xseqMS <- seq(from = xmin, to = xmaxminus, length.out = 1e4)
fxMS <- dMSBPL(x = xseqMS, b = b, C = C, xmin = xmin, xmax = xmax, w = w, v = v)
lines(xseqMS, fxMS, col = mycolors[2], lwd = 2, lty = "dashed") #density curve for bounded power law
legend("topright", bty = "n", lwd = 2, lty = c("dashed"), col = c(mycolors[2]), legend = c("calculated"))

plot(log(xseq), log(fx), type = "l", lwd = 2, lty = "solid", col = mycolors[1], xlab = expression(log(italic(x))), ylab = "log density")
lines(log(xseqMS), log(fxMS), col = mycolors[2], lwd = 2, lty = "dashed") #density curve for bounded power law
legend("topright", bty = "n", lwd = 2, lty = c("solid", "dashed"), col = c(mycolors[1], mycolors[2]), legend = c("bounded power law", "minus-sampled bounded power law"))

plotnegloglik(x = x[isinframe], w = w, v = v)
abline(v = b, lty = "solid", col = mycolors[1])
bML <- estimatebMSBPL(x = x[isinframe], w = w, v = v)
abline(v = bML$minimum, lty = "dashed", col = mycolors[2])
legend("topright", lty = c("solid", "dashed"), legend = c("true", "ML"), col = mycolors[1:2], bty = "n")


