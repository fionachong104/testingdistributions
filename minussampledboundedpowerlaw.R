rm(list = ls())
library(RColorBrewer)

#normalization constant C from Equation 2 in Edwards et al
#arguments:
#min mass xmin
#max mass xmax
#exponent b of bounded power law distribution
#Value: normalization constant C as in Equation 2 of Edwards
getC <- function(xmin, xmax, b){
  if(b == -1){
    return(1 / (log(xmax) - log(xmin)))
  } else {
    return((b + 1) / (xmax ^ (b + 1) - xmin ^(b + 1)))
  }
}

#bounded power law density
#Arguments:
#x: value at which to get density
#b: power law coefficient
#C: normalization constant
#xmin, xmax: min and max for bounded power law
#Value: density f(x) at x
dboundedpowerlaw <- function(x, b, C, xmin, xmax){
  fx <- C * x^b
  fx[x < xmin] <- 0
  fx[x > xmax] <- 0
  return(fx)
}

#calculate pi(r), the probability we minus-sample a circle of radius r
#Arguments:
#w, v: rectangular window dimensions
#r: circle radius
#Value: probability a circle of radius r is minus-sampled
getpir <- function(w, v, r){
  return((w - 2 * r) * (v - 2 * r) / (w * v)) #area in which we can minus-sample a circle of radius r / window area
}

#density for minus-sampled bounded power law
#Arguments:
#x: value at which to get density
#b: power law coefficient
#C: normalization constant
#xmin, xmax: min and max for bounded power law
#w, v: width and height of sampling window (ASSUMES v <= w)
#Value: density f(x) at x
#NOTE: HAVE NOT YET DEALT WITH SPECIAL CASES b = -2, b = -3/2, b = -1
#NOTE: NOT SURE ABOUT JACOBIAN YET
dMSBPL <- function(x, b, C, xmin, xmax, w, v){
  Jacobian <- 1 / (2 * sqrt(pi * x))
  #numerator <- C / (w * v) * (4 / pi * x ^ (b + 1) - 2 * (w + v) / sqrt(pi) * x ^ (b + 0.5) + w * v * x ^ b) #includes factor C / (w * v) that cancels with denominator
  numerator <- 4 / pi * x ^ (b + 1) - 2 * (w + v) / sqrt(pi) * x ^ (b + 0.5) + w * v * x ^ b
  #numerator <- numerator * Jacobian #change of variables needed?
  
  n2 <- dboundedpowerlaw(x = x, b = b, C = C, xmin = xmin, xmax = xmax) * getpir(w = w, v = v, r = sqrt(x / pi))
  #n2 <- n2 * Jacobian #CHECK: JACOBIAN NEEDED HERE?

  xmw <- pi / 4 * v^2 #largest circle we can fit in window
  #dplus <- C / (w * v) * (w * v / (b + 1) * xmw ^ (b + 1) - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * xmw ^ (b + 1.5) + 4 / (pi * (b + 2)) * xmw ^ (b + 2)) #includes factor C / (w * v) that cancels with numerator
  #dminus <- C / (w * v) * (w * v / (b + 1) * xmin ^ (b + 1) - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * xmin ^ (b + 1.5) + 4 / (pi * (b + 2)) * xmin ^ (b + 2)) #includes factor C / (w * v) that cancels with numerator
  if(b == -2){#special case for denominator
    dplus <- w * v / (b + 1) * xmw ^ (b + 1) - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * xmw ^ (b + 1.5) + 4 / (pi * xmw)
    dminus <- w * v / (b + 1) * xmin ^ (b + 1) - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * xmin ^ (b + 1.5) + 4 / (pi * xmin)
    print("b = -2")
  } else if (b == - 1.5){#special case for denominator
    dplus <- w * v / (b + 1) * xmw ^ (b + 1) - 2 * (w + v) / (sqrt(pi) * xmw) + 4 / (pi * (b + 2)) * xmw ^ (b + 2)
    dminus <- w * v / (b + 1) * xmin ^ (b + 1) - 2 * (w + v) / (sqrt(pi) * xmin) + 4 / (pi * (b + 2)) * xmin ^ (b + 2)
    print("b = -1.5")
    print(c(dplus, dminus))
  } else if (b == -1){#special case for denominator
    dplus <- w * v / xmw - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * xmw ^ (b + 1.5) + 4 / (pi * (b + 2)) * xmw ^ (b + 2)
    dminus <- w * v / xmin - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * xmin ^ (b + 1.5) + 4 / (pi * (b + 2)) * xmin ^ (b + 2)
    print("b = -1")
  } else {
    dplus <- w * v / (b + 1) * xmw ^ (b + 1) - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * xmw ^ (b + 1.5) + 4 / (pi * (b + 2)) * xmw ^ (b + 2)
    dminus <- w * v / (b + 1) * xmin ^ (b + 1) - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * xmin ^ (b + 1.5) + 4 / (pi * (b + 2)) * xmin ^ (b + 2)
    print("b: general case")
    print(c(dplus, dminus))
    
  }
  denominator <- dplus - dminus
  d2 <- dMSBPLintegral(b = b, C = C, w = w, v = v, xmin = xmin, xmax = xmax, xmaxminus = xmw)

  fx1 <- numerator / denominator
  fx2 <- n2 / d2$value #USING NUMERICAL INTEGRAL FOR NOW: we can always make this form work for other distributions
  
  print("overall check:")
  print(head(cbind(fx1, fx2)))
  
  fx2[x < xmin] <- 0
  fx2[x > xmax] <- 0
  return(fx2)
}

#numerical check on denominator for minus-sampled bounded power law
#Arguments:
#b: power law coefficient
#C: normalization constant
#w, v: width and height of sampling window (ASSUMES v < w)
#xmin, xmax: min and max for bounded power law
#xmaxminus: max area we can fit in window
#Value: numerical integral of what I think is numerator of minus-sampled bounded power law
dMSBPLintegral <- function(b, C, w, v, xmin, xmax, xmaxminus){
  myfunction <- function(x){
      f <- dboundedpowerlaw(x = x, b = b, C = C, xmin = xmin, xmax = xmax) * getpir(w = w, v = v, r = sqrt(x / pi))
      return(f)
    }
  return(integrate(f = myfunction, lower = xmin, upper = xmaxminus))
}

#inverse cdf for bounded power law, b not equal to -1
#Arguments:
#u in [0, 1]
#xmin and xmax: min and max values of x
#Value: inverse cdf evaluated at u
FXinv <- function(u, b, xmin, xmax){
  (u * xmax ^ (b + 1) + (1 - u) * xmin ^ (b + 1)) ^ (1 / (b + 1))  
}

#simulate from bounded power law distribution
#Arguments:
#n: number of values to simulate
#b: exponent for distribution
#xmin: min size
#xmax: max size
#Value:
#vector of n observations drawn from the specified bounded power law distribution
simulateboundedpowerlaw <- function(n, b, xmin, xmax){
  u <- runif(n = n, min = 0, max = 1)
  x <- FXinv(u = u, b = b, xmin = xmin, xmax = xmax)
}

#radius of circle given area
#Argument: area x
#Value: radius r
radiusfromarea <- function(x){
  return(sqrt(x / pi))
}

#does circle of radius r centred at (alpha, beta) fall entirely within a rectangle of side lengths w, v, lower left corner at origin?
#Arguments:
#alpha, beta: coordinates of circle centre
#r: radius of circle
#w, v: width and height of rectangle
#Value:
#logical: circle entirely in rectangle?
inframe <- function(alpha, beta, r, w, v){
  return(alpha >= r & alpha <= (w - r) & beta >= r & beta <= (v - r))
}

#draw sampling window
#Arguments: w, v width and height of rectangular window
#Value: rectangle representing the window
plotframe <- function(w, v){
  polygon(x = c(0, w, w, 0), y = c(0, 0, v, v))
}

#draw circle, filled if entirely in window
#Arguments:
#alpha, beta: centre of circle
#r: radius of circle
#isinframe (logical): is circle entirely in window?
#Value: draw the circle
plotcircle <- function(alpha, beta, r, isinframe){
  theta <- seq(from = 0, to = 2 * pi, length.out = 1e3)
  xc <- r * cos(theta)
  yc <- r * sin(theta)
  if(isinframe){
    plotcol <- adjustcolor("black", 0.1)
    edgecol <- plotcol
  } else {
    plotcol <- NA
    edgecol <- adjustcolor("black", 0.1)
  }
  polygon(xc + alpha, yc + beta, col = plotcol, border = edgecol)
}

#dimensions of sampling window
w <- 5
v <- 2 #ASSUME v <= w

#load data
#oneyeardf <- read.csv("oneyeardf.csv")
#oneyeardf <- oneyeardf[is.na(oneyeardf$ROI.LabelCode), ] # subset by removing corals that are out of frame
#axisscores <- read.csv("axisscores.csv")
#axisscores <- axisscores[order(axisscores$PC1), ] # order by PC1 from low to high

#sample from the sites; loop to repeat for all sites
#subsample from the minus sampled distribution


#parameters of bounded power law 
# need to make a function that estimates the parameters for each site
xmin <- 1e-1 #min area
xmax <- pi * 3^2 #max area
xmaxminus <- pi / 4 * v^2 #max area of circle that can fit in window
b <- -1.5 #power law exponent
C <- getC(xmin = xmin, xmax = xmax, b = b) #normalization constant
n <- 1e5 #number of circles
nplot <- 100 #number to plot
x <- simulateboundedpowerlaw(n = n, b = b, xmin = xmin, xmax = xmax) #simulated areas

r <- radiusfromarea(x) #colony radii (assume circles)
mycolors <- brewer.pal(3, "Dark2")

alpha <- runif(n = n, min = 0, max = w) #(alpha, beta) uniform random points in sampling window
beta <- runif(n = n, min = 0, max = v)

par(mfrow = c(2, 2))
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
