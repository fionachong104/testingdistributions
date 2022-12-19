rm(list = ls())
library(dplyr)
library(ggplot2)
library(svglite)

# ======================================================================================================
#
# Functions for "Fishing and habitat condition deferentially affect size spectra slopes of coral reef
# fishes". Size spectra functions were adopted from Edwards et al. (2017) --
# https://github.com/andrew-edwards/fitting-size-spectra
#
# ======================================================================================================

# Statistical functions:
# set.params - create list with parameters for running MLE
# negLL.PLB - negative log-likelihood function (function by Edwards et al. 2017)
# pPLB - bounded power-law distribution function (function by Edwards et al. 2017)
# mle_b - # Use analytical value of MLE b for PL model
# slope.conf.int - calculate confidence intervals of estimated slope

# Plotting functions:
# logTicks - add axes and tick marks to a log-log plot to represent
# trophic.plots - plot fitted slope for different trophic positions
# slope_regAndFunc - calculate slopes for each region and trophic group
# slope_reg - calculate slopes for each region

set.params = function(Area){
  # Creates a list of parameters for input into the negll.PLB and pPLB functions
  Area <- Area
  log.Area <- log(Area)
  sum.log.Area <- sum(log.Area)
  min.Area <- min(Area)
  max.Area <- max(Area)
  out.list <- (list(Area, log.Area, sum.log.Area, min.Area, max.Area))
  names(out.list) <- c("Area", "log.Area", "sum.log.Area", "min.Area", "max.Area")
  return(out.list)
}

negLL.PLB = function(b, x, n, xmin, xmax, sumlogx)
  {
    if(xmin <= 0 | xmin >= xmax) stop("Parameters out of bounds in negLL.PLB")
    if(b != -1)
      { neglogLL = -n * log( ( b + 1) / (xmax^(b + 1) - xmin^(b + 1)) ) -
            b * sumlogx
      } else
      { neglogLL = n * log( log(xmax) - log(xmin) ) + sumlogx
      }
    return(neglogLL)
}

pPLB = function(x = 10, b = -2, xmin = 1, xmax = 100)    
  {
    if(xmin <= 0 | xmin >= xmax) stop("Parameters out of bounds in pPLB")
    y = 0 * x     # so have zeros where x < xmin
    y[x > xmax] = 1  # 1 for x > xmax
    if(b != -1)
        {  xmintobplus1 = xmin^(b+1)
           denom = xmax^(b+1) - xmintobplus1
           y[x >= xmin & x <= xmax] =
               ( x[x >= xmin & x <= xmax]^(b + 1) - xmintobplus1 ) / denom
        } else
        {  logxmin = log(xmin)
           denom = log(xmax) - logxmin
           y[x >= xmin & x <= xmax] =
               ( log( x[x >= xmin & x <= xmax] ) - logxmin ) / denom
        }
    return(y)                              
  }

mle_b = function(region, x, log_x, sum_log_x, x_min, x_max){
    # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2017)
    # as a starting point for nlm for MLE of b for PLB model. Code adopted from
    # Edwards et al. (2017).

    PL.bMLE = 1/( log(min(x)) - sum_log_x/length(x)) - 1
    
    PLB.minLL =  nlm(negLL.PLB, p=PL.bMLE, x=x, n=length(x),
        xmin=x_min, xmax=x_max, sumlogx=sum_log_x, hessian = TRUE) #, print.level=2 )
    
    PLB.bMLE = PLB.minLL$estimate
    
    PLB.return = list(PLB.bMLE, PLB.minLL)
    
    return(PLB.return)
}

slope.conf.int = function(PLB.bMLE.b, PLB.minLL.b, input){
  # Calculate confidence intervals. Code adopted from Edwards et al. (2017)
  bvec = seq(PLB.bMLE.b - 0.5, PLB.bMLE.b + 0.5, 0.00001) 
  PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
  for(i in 1:length(bvec)){
    PLB.LLvals[i] = negLL.PLB(bvec[i], x=input$Area, n=length(input$Area), xmin=input$min.Area,
        xmax=input$max.Area, sumlogx=input$sum.log.Area)   
  }
  critVal = PLB.minLL.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
  bIn95 = bvec[ PLB.LLvals < critVal ]
  return(c(min(bIn95), max(bIn95)))
}

logTicks = function(xLim, yLim = NULL, tclSmall = -0.2, xLabelSmall = NULL,
                    yLabelSmall = NULL, xLabelBig = NULL, mgpVal=c(1.6,0.5,0)){
    # Modified version of function by Edwards et al. (2017)
    ll = 1:9
    log10ll = log10(ll)
    # box()                                                                       # PC removed
    # x axis
    if(!is.null(xLim))               # if NULL then ignore
      {  
      # Do enough tick marks to encompass axes:
      xEncompassLog = c(floor(log10(xLim[1])), ceiling(log10(xLim[2])))
      xBig = 10^c(xEncompassLog[1]:xEncompassLog[2])
      # Big unlabelled, always want these:
      # axis(1, at= xBig, labels = rep("", length(xBig)), mgp = mgpVal)           # PC removed
      # Big labelled:
      if(is.null(xLabelBig)) { xLabelBig = xBig }
      axis(1, at= xLabelBig, labels = xLabelBig, mgp = mgpVal)
      # axis(1, at=c(1, 10, 100), labels = c(1, 10, 100), mgp=c(1.7,0.7,0))
      # Small unlabelled:
      # axis(1, xBig %x% ll, labels=rep("", length(xBig %x% ll)), tcl=tclSmall)
      # Small labelled:
      if(!is.null(xLabelSmall))
          {
          axis(1, at=xLabelSmall, labels=xLabelSmall, mgp=mgpVal, tcl=tclSmall)
          }
      }
    # Repeat for y axis:
    if(!is.null(yLim))               # if NULL then ignore
      {  
      # Do enough tick marks to encompass axes:
      yEncompassLog = c(floor(log10(yLim[1])), ceiling(log10(yLim[2])))
      # yBig = 10^c(yEncompassLog[1]:yEncompassLog[2])                            # PC removed
      yBig = c(0,1,10,100,1000,10000)
      # Big labelled:
      axis(2, at= yBig, labels = yBig, mgp = mgpVal)
      # Small unlabelled:
      # axis(2, yBig %x% ll, labels=rep("", length(yBig %x% ll)), tcl=tclSmall)   # PC removed
      # Small labelled:
      if(!is.null(yLabelSmall))
          {
          axis(2, at=yLabelSmall, labels=yLabelSmall, mgp=mgpVal, tcl=tclSmall)
          }
      }     
}

trophic.plots = function(PLB.return, PLB.bMLE.b, PLB.minLL.b, df.input, mgpVals, troph_id, panel){
  # plot and find 95% confidence intervals for MLE method.
  PLB.minNegLL.b <- PLB.minLL.b$minimum
  x <- df.input$Area  
  plot(sort(df.input$Area, decreasing=TRUE), 1:length(df.input$Area), log="xy",
    xlab=expression(paste("Area, ", italic(x), " (log[cm])")),
    # xlab = "log(kg)",
    ylab = expression(paste("Number of corals", " ">=" ", italic("x"))), mgp=mgpVals,
    # ylab = "  ",
    xlim = c(df.input$min.Area, df.input$max.Area), ylim = c(1, 10000), axes=FALSE)
  logTicks(xLim, yLim, xLabelBig = c(0, 1, 5, 10)) # Tick marks.
  x.PLB = seq(min(df.input$Area), max(df.input$Area), length=1000) # x values to plot PLB. Note
                                                                         # that these encompass the data, and are not based
                                                                         # on the binning (in MEE Figure 6 the line starts as
                                                                         # min(x), not the first bin.
  y.PLB = (1 - pPLB(x = x.PLB, b = PLB.bMLE.b, xmin = min(x.PLB),
    xmax = max(x.PLB))) * length(df.input$Area)
  lines(x.PLB, y.PLB, col="red", lwd=2)
  text(x = 0.1, y = 15, labels = troph_id, cex = 1.1, pos = 1, col = "black")
  spectra.text <- as.character(round(PLB.bMLE.b, 2))
  text(x=0.1, y=5, labels = bquote(paste(italic("b = "), .(spectra.text))), 
    cex=1.1, pos=1, col="black")
  mtext(panel, side = 3, cex = 1.4, adj = -0.15)
  # Values of b to test to obtain confidence interval. For the real movement data
  # sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
  # symmetric interval here.
  bvec = seq(PLB.bMLE.b - 0.5, PLB.bMLE.b + 0.5, 0.00001) 
  PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
  for(i in 1:length(bvec)){
    PLB.LLvals[i] = negLL.PLB(bvec[i], x=df.input$Area, n=length(df.input$Area), xmin=df.input$min.Area,
      xmax=df.input$max.Area, sumlogx=df.input$sum.log.Area)   
  }
  critVal = PLB.minNegLL.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
  bIn95 = bvec[ PLB.LLvals < critVal ]
  # To add just the curves at the limits of the 95% confidence interval of b:
  for(i in c(1, length(bIn95))){
    lines(x.PLB, (1 - pPLB(x = x.PLB, b = bIn95[i], xmin = min(x.PLB),
      xmax = max(x.PLB))) * length(df.input$Area), col="red", lty=2)
  }
}

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
  fx <- C * x ^ b
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

#density for minus-sampled bounded power law (in form we can supply to optimizer)
#Arguments:
#b: power law coefficient
#x: value at which to get density
#xmin, min for bounded power law
#w, v: width and height of sampling window (ASSUMES v <= w and xmax greater than largest observable object)
#Value: density f(x) at x
dMSBPLanalytic <- function(b, x, xmin, w, v){
  numerator <- 4 / pi * x ^ (b + 1) - 2 * (w + v) / sqrt(pi) * x ^ (b + 0.5) + w * v * x ^ b
  xmw <- pi / 4 * v ^ 2 #largest circle we can fit in window
  if(b == -2){#special case for denominator
    dplus <- w * v / (b + 1) * xmw ^ (b + 1) - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * xmw ^ (b + 1.5) + 4 / pi * log(xmw)
    dminus <- w * v / (b + 1) * xmin ^ (b + 1) - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * xmin ^ (b + 1.5) + 4 / pi * log(xmin)
  } else if (b == - 1.5){#special case for denominator
    dplus <- w * v / (b + 1) * xmw ^ (b + 1) - 2 * (w + v) / sqrt(pi) * log(xmw) + 4 / (pi * (b + 2)) * xmw ^ (b + 2)
    dminus <- w * v / (b + 1) * xmin ^ (b + 1) - 2 * (w + v) / sqrt(pi) * log(xmin) + 4 / (pi * (b + 2)) * xmin ^ (b + 2)
  } else if (b == -1){#special case for denominator
    dplus <- w * v * log(xmw) - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * xmw ^ (b + 1.5) + 4 / (pi * (b + 2)) * xmw ^ (b + 2)
    dminus <- w * v * log(xmin) - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * xmin ^ (b + 1.5) + 4 / (pi * (b + 2)) * xmin ^ (b + 2)
  } else {
    dplus <- w * v / (b + 1) * xmw ^ (b + 1) - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * xmw ^ (b + 1.5) + 4 / (pi * (b + 2)) * xmw ^ (b + 2)
    dminus <- w * v / (b + 1) * xmin ^ (b + 1) - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * xmin ^ (b + 1.5) + 4 / (pi * (b + 2)) * xmin ^ (b + 2)
  }
  denominator <- dplus - dminus
  fx <- numerator / denominator
  fx[x < xmin] <- 0
  fx[x > xmw] <- 0
  return(fx)
}

#negative log likelihood for minus-sampled bounded power law (in form we can supply to optimizer)
#Arguments:
#b: power law coefficient
#x: vector of sizes
#xmin, min for bounded power law
#w, v: width and height of sampling window (ASSUMES v <= w and xmax greater than largest observable object)
#Value: negative log likelihood for observations x, with parameter b
negloglikMSBPL <- function(b, x, w, v){
  xmin <- min(x) #this is maximum likelihood estimate
  logfx <- sum(log(dMSBPLanalytic(b = b, x = x, xmin = xmin, w = w, v = v)))
  return(-logfx)
}

#plot negative log likelihood against b for minus-sampled bounded power law
#Arguments:
#x: vector of sizes
#w: width of window
#v: height of window
#Value: plot of negative log likelihood against b
plotnegloglik <- function(x, w, v){
  b <- seq(from = -3, to = -1e-6, length.out = 100)
  negllvec <- numeric(3)
  for(i in 1:100){
    negllvec[i] <- negloglikMSBPL(b = b[i], x = x, w = w, v = v)
  }
  plot(b, negllvec, type = "l", xlab = expression(italic(b)), ylab = expression(-log(italic(f(x)))))  
}

#maximum likelihood estimate of b for minus-sampled bounded power law
#Arguments:
#x: vector of sizes
#w: width of window
#v: height of window
#Value: object returned by optimize(), for which $minimum is the ML estimate and $objective is the negative log likelihood
estimatebMSBPL <- function(x, w, v){
  return(optimize(f = negloglikMSBPL, interval = c(-3, 0), x = x, w = w, v = v))
}

#density for minus-sampled bounded power law
#Arguments:
#x: value at which to get density
#b: power law coefficient
#C: normalization constant
#xmin, xmax: min and max for bounded power law
#w, v: width and height of sampling window (ASSUMES v <= w)
#Value: density f(x) at x
#NOTE: NOT SURE ABOUT JACOBIAN YET
dMSBPL <- function(x, b, C, xmin, xmax, w, v){
  Jacobian <- 1 / (2 * sqrt(pi * x))
  #numerator <- C / (w * v) * (4 / pi * x ^ (b + 1) - 2 * (w + v) / sqrt(pi) * x ^ (b + 0.5) + w * v * x ^ b) #includes factor C / (w * v) that cancels with denominator
  numerator <- 4 / pi * x ^ (b + 1) - 2 * (w + v) / sqrt(pi) * x ^ (b + 0.5) + w * v * x ^ b
  #numerator <- numerator * Jacobian #change of variables needed?
  
  n2 <- dboundedpowerlaw(x = x, b = b, C = C, xmin = xmin, xmax = xmax) * getpir(w = w, v = v, r = sqrt(x / pi))
  #n2 <- n2 * Jacobian #CHECK: JACOBIAN NEEDED HERE?
  
  xmw <- pi / 4 * v ^ 2 #largest circle we can fit in window
  if(b == -2){#special case for denominator
    dplus <- w * v / (b + 1) * xmw ^ (b + 1) - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * xmw ^ (b + 1.5) + 4 / pi * log(xmw)
    dminus <- w * v / (b + 1) * xmin ^ (b + 1) - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * xmin ^ (b + 1.5) + 4 / pi * log(xmin)
  } else if (b == - 1.5){#special case for denominator
    dplus <- w * v / (b + 1) * xmw ^ (b + 1) - 2 * (w + v) / sqrt(pi) * log(xmw) + 4 / (pi * (b + 2)) * xmw ^ (b + 2)
    dminus <- w * v / (b + 1) * xmin ^ (b + 1) - 2 * (w + v) / sqrt(pi) * log(xmin) + 4 / (pi * (b + 2)) * xmin ^ (b + 2)
  } else if (b == -1){#special case for denominator
    dplus <- w * v * log(xmw) - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * xmw ^ (b + 1.5) + 4 / (pi * (b + 2)) * xmw ^ (b + 2)
    dminus <- w * v * log(xmin) - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * xmin ^ (b + 1.5) + 4 / (pi * (b + 2)) * xmin ^ (b + 2)
  } else {
    dplus <- w * v / (b + 1) * xmw ^ (b + 1) - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * xmw ^ (b + 1.5) + 4 / (pi * (b + 2)) * xmw ^ (b + 2)
    dminus <- w * v / (b + 1) * xmin ^ (b + 1) - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * xmin ^ (b + 1.5) + 4 / (pi * (b + 2)) * xmin ^ (b + 2)
  }
  denominator <- dplus - dminus
  d2 <- dMSBPLintegral(b = b, C = C, w = w, v = v, xmin = xmin, xmax = xmax, xmaxminus = xmw)
  
  fx1 <- numerator / denominator
  fx2 <- n2 / d2$value #USING NUMERICAL INTEGRAL FOR NOW: we can always make this form work for other distributions
  
  print("overall check:")
  print(head(cbind(fx1, fx2)))
  
  fx2[x < xmin] <- 0
  fx2[x > xmw] <- 0
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

#inverse cdf for bounded power law
#Arguments:
#u in [0, 1]
#xmin and xmax: min and max values of x
#Value: inverse cdf evaluated at u
FXinv <- function(u, b, xmin, xmax){
  if(b == -1){
    x <- xmax ^ u * xmin ^ (1 - u)
  } else {
    x <- (u * xmax ^ (b + 1) + (1 - u) * xmin ^ (b + 1)) ^ (1 / (b + 1))  
  }
  return(x)
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
# w <- 5
# v <- 2 #ASSUME v <= w

# images dimensions from CANON camera (w x v = 3648 x 2736)
# 10 cm = roughly 350 pixels; 1 m = 3500
w <- 3648/3500
v <- 2736/3500



#load data

oneyeardf <- read.csv("oneyeardf.csv",row.names=1)
oneyeardf <- oneyeardf[is.na(oneyeardf$ROI.LabelCode), ] # to remove corals that were out of frame 
axisscores <- read.csv("axisscores.csv", row.names=1)
axisscores <- axisscores[order(axisscores$PC1), ]

# #Set MLE parameters ===================================================================================
#   # Areas for each site
# 
# BR <- oneyeardf %>% filter(Site== "Black Rock")
# CI <- oneyeardf %>% filter(Site== "Cook Island")
# FR <- oneyeardf %>% filter(Site== "Flat Rock")
# FLI <- oneyeardf %>% filter(Site== "Flinders")
# IGS <- oneyeardf %>% filter(Site== "Inner Gneering Shoals")
# Hen <- oneyeardf %>% filter(Site== "Henderson Rock")
# JRFT <- oneyeardf %>% filter(Site== "Julian Rock False Trench")
# JRN <- oneyeardf %>% filter(Site== "Julian Rock Nursery")
# LEI <- oneyeardf %>% filter(Site== "Lady Elliot Island")
# LMI <- oneyeardf %>% filter(Site== "Lady Musgrave Island")
# LL <- oneyeardf %>% filter(Site== "Libbies Lair")
# Mud <- oneyeardf %>% filter(Site== "Mudjimba")
# NR <- oneyeardf %>% filter(Site== "North Rock")
# NSI <- oneyeardf %>% filter(Site== "North Solitary Island")
# NWSI <- oneyeardf %>% filter(Site== "North West Solitary Island")
# SSI <- oneyeardf %>% filter(Site== "South Solitary Island")
# SWSI <- oneyeardf %>% filter(Site== "South West Solitary Island")
# Ten <- oneyeardf %>% filter(Site== "Tenements")
# WR <- oneyeardf %>% filter(Site== "Wolf Rock")
# Wool <- oneyeardf %>% filter(Site== "Woolgoolga Reef")
# 
# BR.input <- set.params(BR$Area)
# CI.input <- set.params(CI$Area)
# FR.input <- set.params(FR$Area)
# FLI.input <- set.params(FLI$Area)
# IGS.input <- set.params(IGS$Area)
# Hen.input <- set.params(Hen$Area)
# JRFT.input <- set.params(JRFT$Area)
# JRN.input <- set.params(JRN$Area)
# LEI.input <- set.params(LEI$Area)
# LMI.input <- set.params(LMI$Area)
# LL.input <- set.params(LL$Area)
# Mud.input <- set.params(Mud$Area)
# NR.input <- set.params(NR$Area)
# NSI.input <- set.params(NSI$Area)
# NWSI.input <- set.params(NWSI$Area)
# SSI.input <- set.params(SSI$Area)
# SWSI.input <- set.params(SWSI$Area)
# Ten.input <- set.params(Ten$Area)
# WR.input <- set.params(WR$Area)
# Wool.input <- set.params(Wool$Area)
# 
# 
# mgpVals <- c(1.6,0.5,0) # mgp values 2.0, 0.5, 0
# xLim <- 10^par("usr")[1:2]
# yLim <- 10^par("usr")[3:4]
# 
# # MLE Henderson Rock (Hen) Area ===============================================================================
# # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# # as a starting point for nlm for MLE of b for PLB model.
# PLB.return.Hen <- mle_b(Site== "Henderson Rock", x=Hen.input$Area, log_x=Hen.input$log.Area, sum_log_x=Hen.input$sum.log.Area,
#                        x_min=Hen.input$min.Area, x_max=Hen.input$max.Area)
# PLB.bMLE.Hen.b <- PLB.return.Hen[[1]] 
# PLB.minLL.Hen.b <- PLB.return.Hen[[2]]
# 
# # plot and find 95% confidence intervals for MLE method.
# PLB.minNegLL.Hen.b <- PLB.minLL.Hen.b$minimum
# x <- Hen.input$Area
# Henx.PLB = seq(min(Hen.input$Area), max(Hen.input$Area), length=1000) # x values to plot PLB. Note
# # that these encompass the data, and are not based
# # on the binning (in MEE Figure 6 the line starts as
# # min(x), not the first bin.
# Heny.PLB = (1 - pPLB(x = Henx.PLB, b = PLB.bMLE.Hen.b, xmin = min(Henx.PLB),
#                     xmax = max(Henx.PLB))) * length(Hen.input$Area)
# spectra.text <- as.character(round(PLB.bMLE.Hen.b, 2))
# Henb_plot <- ggplot() +
#   geom_point(aes(x = (sort(Hen.input$Area, decreasing=TRUE)), y = (1:length(Hen.input$Area))),
#              color = "cadetblue", size = 2, alpha = 0.3) +
#   xlab(expression(paste("Colony area, ", italic("x"), " (cm^2)"))) +
#   ylab(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    "))) +
#   scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
#                      limits = c(0.25, max(length(Hen.input$Area),
#                                           length(BR.input$Area),
#                                           length(CI.input$Area),
#                                           length(FR.input$Area),
#                                           length(FLI.input$Area),
#                                           length(IGS.input$Area),
#                                           length(JRFT.input$Area),
#                                           length(JRN.input$Area),
#                                           length(LEI.input$Area),
#                                           length(LMI.input$Area),
#                                           length(LL.input$Area),
#                                           length(Mud.input$Area),
#                                           length(NSI.input$Area),
#                                           length(NWSI.input$Area),
#                                           length(SSI.input$Area),
#                                           length(SWSI.input$Area),
#                                           length(Ten.input$Area),
#                                           length(WR.input$Area),
#                                           length(Wool.input$Area)))) +
#   scale_x_continuous(trans = 'log10', breaks = c(0,1,5,10,100,1000,10000),
#                      limits = c(min(Hen.input$min.Area,
#                                     BR.input$min.Area,
#                                     CI.input$min.Area,
#                                     FR.input$min.Area,
#                                     FLI.input$min.Area,
#                                     IGS.input$min.Area,
#                                     JRFT.input$min.Area,
#                                     JRN.input$min.Area,
#                                     LEI.input$min.Area,
#                                     LMI.input$min.Area,
#                                     LL.input$min.Area,
#                                     Mud.input$min.Area,
#                                     NSI.input$min.Area,
#                                     NWSI.input$min.Area,
#                                     SSI.input$min.Area,
#                                     SWSI.input$min.Area,
#                                     Ten.input$min.Area,
#                                     WR.input$min.Area,
#                                     Wool.input$min.Area), 
#                                 max(Hen.input$max.Area,
#                                     BR.input$max.Area,
#                                     CI.input$max.Area,
#                                     FR.input$max.Area,
#                                     FLI.input$max.Area,
#                                     IGS.input$max.Area,
#                                     JRFT.input$max.Area,
#                                     JRN.input$max.Area,
#                                     LEI.input$max.Area,
#                                     LMI.input$max.Area,
#                                     LL.input$max.Area,
#                                     Mud.input$max.Area,
#                                     NSI.input$max.Area,
#                                     NWSI.input$max.Area,
#                                     SSI.input$max.Area,
#                                     SWSI.input$max.Area,
#                                     Ten.input$max.Area,
#                                     WR.input$max.Area,
#                                     Wool.input$max.Area))) +
#   geom_line(aes(x = Henx.PLB, y = Heny.PLB), col = 'black', lwd = 1) +
#   annotate("text", x=5, y=10, label="Henderson Rock") +
#   annotate("text", x=5, y=3, label = expression(paste(italic("b = "), -0.929))) +
#   theme_classic()
# 
# # Values of b to test to obtain confidence interval. For the real movement data
# # sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# # symmetric interval here.
# bvec = seq(PLB.bMLE.Hen.b - 0.5, PLB.bMLE.Hen.b + 0.5, 0.00001) 
# PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
# for(i in 1:length(bvec)){
#   PLB.LLvals[i] = negLL.PLB(bvec[i], x=Hen.input$Area, n=length(Hen.input$Area), xmin=Hen.input$min.Area,
#                             xmax=Hen.input$max.Area, sumlogx=Hen.input$sum.log.Area)   
# }
# critVal = PLB.minNegLL.Hen.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
# bIn95 = bvec[ PLB.LLvals < critVal ]
# HenbIn95 <- c(min(bIn95), max(bIn95))
# 
# # MLE North Rock (NR) Area ===============================================================================
# # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# # as a starting point for nlm for MLE of b for PLB model.
# PLB.return.NR <- mle_b(Site== "North Rock", x=NR.input$Area, log_x=NR.input$log.Area, sum_log_x=NR.input$sum.log.Area,
#                         x_min=NR.input$min.Area, x_max=NR.input$max.Area)
# PLB.bMLE.NR.b <- PLB.return.NR[[1]] 
# PLB.minLL.NR.b <- PLB.return.NR[[2]]
# 
# # plot and find 95% confidence intervals for MLE method.
# PLB.minNegLL.NR.b <- PLB.minLL.NR.b$minimum
# x <- NR.input$Area
# NRx.PLB = seq(min(NR.input$Area), max(NR.input$Area), length=1000) # x values to plot PLB. Note
# # that these encompass the data, and are not based
# # on the binning (in MEE Figure 6 the line starts as
# # min(x), not the first bin.
# NRy.PLB = (1 - pPLB(x = NRx.PLB, b = PLB.bMLE.NR.b, xmin = min(NRx.PLB),
#                      xmax = max(NRx.PLB))) * length(NR.input$Area)
# spectra.text <- as.character(round(PLB.bMLE.NR.b, 2))
# NRb_plot <- ggplot() +
#   geom_point(aes(x = (sort(NR.input$Area, decreasing=TRUE)), y = (1:length(NR.input$Area))),
#              color = "cadetblue", size = 2, alpha = 0.3) +
#   xlab(expression(paste("Colony area, ", italic("x"), " (cm^2)"))) +
#   ylab(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    "))) +
#   scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
#                      limits = c(0.25, max(length(Hen.input$Area),
#                                           length(NR.input$Area),
#                                           length(BR.input$Area),
#                                           length(CI.input$Area),
#                                           length(FR.input$Area),
#                                           length(FLI.input$Area),
#                                           length(IGS.input$Area),
#                                           length(JRFT.input$Area),
#                                           length(JRN.input$Area),
#                                           length(LEI.input$Area),
#                                           length(LMI.input$Area),
#                                           length(LL.input$Area),
#                                           length(Mud.input$Area),
#                                           length(NSI.input$Area),
#                                           length(NWSI.input$Area),
#                                           length(SSI.input$Area),
#                                           length(SWSI.input$Area),
#                                           length(Ten.input$Area),
#                                           length(WR.input$Area),
#                                           length(Wool.input$Area)))) +
#   scale_x_continuous(trans = 'log10', breaks = c(0,1,5,10,100,1000,10000),
#                      limits = c(min(Hen.input$min.Area,
#                                     NR.input$min.Area,
#                                     BR.input$min.Area,
#                                     CI.input$min.Area,
#                                     FR.input$min.Area,
#                                     FLI.input$min.Area,
#                                     IGS.input$min.Area,
#                                     JRFT.input$min.Area,
#                                     JRN.input$min.Area,
#                                     LEI.input$min.Area,
#                                     LMI.input$min.Area,
#                                     LL.input$min.Area,
#                                     Mud.input$min.Area,
#                                     NSI.input$min.Area,
#                                     NWSI.input$min.Area,
#                                     SSI.input$min.Area,
#                                     SWSI.input$min.Area,
#                                     Ten.input$min.Area,
#                                     WR.input$min.Area,
#                                     Wool.input$min.Area), 
#                                 max(Hen.input$max.Area,
#                                     NR.input$max.Area,
#                                     BR.input$max.Area,
#                                     CI.input$max.Area,
#                                     FR.input$max.Area,
#                                     FLI.input$max.Area,
#                                     IGS.input$max.Area,
#                                     JRFT.input$max.Area,
#                                     JRN.input$max.Area,
#                                     LEI.input$max.Area,
#                                     LMI.input$max.Area,
#                                     LL.input$max.Area,
#                                     Mud.input$max.Area,
#                                     NSI.input$max.Area,
#                                     NWSI.input$max.Area,
#                                     SSI.input$max.Area,
#                                     SWSI.input$max.Area,
#                                     Ten.input$max.Area,
#                                     WR.input$max.Area,
#                                     Wool.input$max.Area))) +
#   geom_line(aes(x = NRx.PLB, y = NRy.PLB), col = 'black', lwd = 1) +
#   annotate("text", x=5, y=10, label="North Rock") +
#   annotate("text", x=5, y=3, label = expression(paste(italic("b = "), -0.100))) +
#   theme_classic()
# 
# # Values of b to test to obtain confidence interval. For the real movement data
# # sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# # symmetric interval here.
# bvec = seq(PLB.bMLE.NR.b - 0.5, PLB.bMLE.NR.b + 0.5, 0.00001) 
# PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
# for(i in 1:length(bvec)){
#   PLB.LLvals[i] = negLL.PLB(bvec[i], x=NR.input$Area, n=length(NR.input$Area), xmin=NR.input$min.Area,
#                             xmax=NR.input$max.Area, sumlogx=NR.input$sum.log.Area)   
# }
# critVal = PLB.minNegLL.NR.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
# bIn95 = bvec[ PLB.LLvals < critVal ]
# NRbIn95 <- c(min(bIn95), max(bIn95))
# 
# # MLE Black Rock (BR) Area ===============================================================================
# # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# # as a starting point for nlm for MLE of b for PLB model.
# PLB.return.BR <- mle_b(Site== "Black Rock", x=BR.input$Area, log_x=BR.input$log.Area, sum_log_x=BR.input$sum.log.Area,
#                        x_min=BR.input$min.Area, x_max=BR.input$max.Area)
# PLB.bMLE.BR.b <- PLB.return.BR[[1]] 
# PLB.minLL.BR.b <- PLB.return.BR[[2]]
# 
# # plot and find 95% confidence intervals for MLE method.
# PLB.minNegLL.BR.b <- PLB.minLL.BR.b$minimum
# x <- BR.input$Area
# BRx.PLB = seq(min(BR.input$Area), max(BR.input$Area), length=1000) # x values to plot PLB. Note
# # that these encompass the data, and are not based
# # on the binning (in MEE Figure 6 the line starts as
# # min(x), not the first bin.
# BRy.PLB = (1 - pPLB(x = BRx.PLB, b = PLB.bMLE.BR.b, xmin = min(BRx.PLB),
#                     xmax = max(BRx.PLB))) * length(BR.input$Area)
# spectra.text <- as.character(round(PLB.bMLE.BR.b, 2))
# BRb_plot <- ggplot() +
#   geom_point(aes(x = (sort(BR.input$Area, decreasing=TRUE)), y = (1:length(BR.input$Area))),
#              color = "cadetblue", size = 2, alpha = 0.3) +
#   xlab(expression(paste("Colony area, ", italic("x"), " (cm^2)"))) +
#   ylab(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    "))) +
#   scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
#                      limits = c(0.25, max(length(BR.input$Area),
#                                           length(CI.input$Area),
#                                           length(FR.input$Area),
#                                           length(FLI.input$Area),
#                                           length(IGS.input$Area),
#                                           length(JRFT.input$Area),
#                                           length(JRN.input$Area),
#                                           length(LEI.input$Area),
#                                           length(LMI.input$Area),
#                                           length(LL.input$Area),
#                                           length(Mud.input$Area),
#                                           length(NSI.input$Area),
#                                           length(NWSI.input$Area),
#                                           length(SSI.input$Area),
#                                           length(SWSI.input$Area),
#                                           length(Ten.input$Area),
#                                           length(WR.input$Area),
#                                           length(Wool.input$Area)))) +
#   scale_x_continuous(trans = 'log10', breaks = c(0,1,5,10,100,1000,10000),
#                      limits = c(min(BR.input$min.Area,
#                                     CI.input$min.Area,
#                                     FR.input$min.Area,
#                                     FLI.input$min.Area,
#                                     IGS.input$min.Area,
#                                     JRFT.input$min.Area,
#                                     JRN.input$min.Area,
#                                     LEI.input$min.Area,
#                                     LMI.input$min.Area,
#                                     LL.input$min.Area,
#                                     Mud.input$min.Area,
#                                     NSI.input$min.Area,
#                                     NWSI.input$min.Area,
#                                     SSI.input$min.Area,
#                                     SWSI.input$min.Area,
#                                     Ten.input$min.Area,
#                                     WR.input$min.Area,
#                                     Wool.input$min.Area), 
#                                 max(BR.input$max.Area,
#                                     CI.input$max.Area,
#                                     FR.input$max.Area,
#                                     FLI.input$max.Area,
#                                     IGS.input$max.Area,
#                                     JRFT.input$max.Area,
#                                     JRN.input$max.Area,
#                                     LEI.input$max.Area,
#                                     LMI.input$max.Area,
#                                     LL.input$max.Area,
#                                     Mud.input$max.Area,
#                                     NSI.input$max.Area,
#                                     NWSI.input$max.Area,
#                                     SSI.input$max.Area,
#                                     SWSI.input$max.Area,
#                                     Ten.input$max.Area,
#                                     WR.input$max.Area,
#                                     Wool.input$max.Area))) +
#   geom_line(aes(x = BRx.PLB, y = BRy.PLB), col = 'black', lwd = 1) +
#   annotate("text", x=5, y=10, label="Black Rock") +
#   annotate("text", x=5, y=3, label = expression(paste(italic("b = "), -0.858))) +
#   theme_classic()
# 
# # Values of b to test to obtain confidence interval. For the real movement data
# # sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# # symmetric interval here.
# bvec = seq(PLB.bMLE.BR.b - 0.5, PLB.bMLE.BR.b + 0.5, 0.00001) 
# PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
# for(i in 1:length(bvec)){
#   PLB.LLvals[i] = negLL.PLB(bvec[i], x=BR.input$Area, n=length(BR.input$Area), xmin=BR.input$min.Area,
#                             xmax=BR.input$max.Area, sumlogx=BR.input$sum.log.Area)   
# }
# critVal = PLB.minNegLL.BR.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
# bIn95 = bvec[ PLB.LLvals < critVal ]
# BRbIn95 <- c(min(bIn95), max(bIn95))
# 
# 
# # MLE Cook Island (CI) Area ===============================================================================
# # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# # as a starting point for nlm for MLE of b for PLB model.
# PLB.return.CI <- mle_b(Site== "Cook Island", x=CI.input$Area, log_x=CI.input$log.Area, sum_log_x=CI.input$sum.log.Area,
#                        x_min=CI.input$min.Area, x_max=CI.input$max.Area)
# PLB.bMLE.CI.b <- PLB.return.CI[[1]] 
# PLB.minLL.CI.b <- PLB.return.CI[[2]]
# 
# # plot and find 95% confidence intervals for MLE method.
# PLB.minNegLL.CI.b <- PLB.minLL.CI.b$minimum
# x <- CI.input$Area
# CIx.PLB = seq(min(CI.input$Area), max(CI.input$Area), length=1000) # x values to plot PLB. Note
# # that these encompass the data, and are not based
# # on the binning (in MEE Figure 6 the line starts as
# # min(x), not the first bin.
# CIy.PLB = (1 - pPLB(x = CIx.PLB, b = PLB.bMLE.CI.b, xmin = min(CIx.PLB),
#                     xmax = max(CIx.PLB))) * length(CI.input$Area)
# spectra.text <- as.character(round(PLB.bMLE.CI.b, 2))
# CIb_plot <- ggplot() +
#   geom_point(aes(x = (sort(CI.input$Area, decreasing=TRUE)), y = (1:length(CI.input$Area))),
#              color = "cadetblue", size = 2, alpha = 0.3) +
#   xlab(expression(paste("Colony area, ", italic("x"), " (cm^2)"))) +
#   ylab(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    "))) +
#   scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
#                      limits = c(0.25, max(length(BR.input$Area),
#                                           length(CI.input$Area),
#                                           length(FR.input$Area),
#                                           length(FLI.input$Area),
#                                           length(IGS.input$Area),
#                                           length(JRFT.input$Area),
#                                           length(JRN.input$Area),
#                                           length(LEI.input$Area),
#                                           length(LMI.input$Area),
#                                           length(LL.input$Area),
#                                           length(Mud.input$Area),
#                                           length(NSI.input$Area),
#                                           length(NWSI.input$Area),
#                                           length(SSI.input$Area),
#                                           length(SWSI.input$Area),
#                                           length(Ten.input$Area),
#                                           length(WR.input$Area),
#                                           length(Wool.input$Area)))) +
#   scale_x_continuous(trans = 'log10', breaks = c(0,1,5,10,100,1000,10000),
#                      limits = c(min(BR.input$min.Area,
#                                     CI.input$min.Area,
#                                     FR.input$min.Area,
#                                     FLI.input$min.Area,
#                                     IGS.input$min.Area,
#                                     JRFT.input$min.Area,
#                                     JRN.input$min.Area,
#                                     LEI.input$min.Area,
#                                     LMI.input$min.Area,
#                                     LL.input$min.Area,
#                                     Mud.input$min.Area,
#                                     NSI.input$min.Area,
#                                     NWSI.input$min.Area,
#                                     SSI.input$min.Area,
#                                     SWSI.input$min.Area,
#                                     Ten.input$min.Area,
#                                     WR.input$min.Area,
#                                     Wool.input$min.Area), 
#                                 max(BR.input$max.Area,
#                                     CI.input$max.Area,
#                                     FR.input$max.Area,
#                                     FLI.input$max.Area,
#                                     IGS.input$max.Area,
#                                     JRFT.input$max.Area,
#                                     JRN.input$max.Area,
#                                     LEI.input$max.Area,
#                                     LMI.input$max.Area,
#                                     LL.input$max.Area,
#                                     Mud.input$max.Area,
#                                     NSI.input$max.Area,
#                                     NWSI.input$max.Area,
#                                     SSI.input$max.Area,
#                                     SWSI.input$max.Area,
#                                     Ten.input$max.Area,
#                                     WR.input$max.Area,
#                                     Wool.input$max.Area))) +
#   geom_line(aes(x = CIx.PLB, y = CIy.PLB), col = 'black', lwd = 1) +
#   annotate("text", x=5, y=10, label="Cook Island") +
#   annotate("text", x=5, y=3, label = expression(paste(italic("b = "), -0.989))) +
#   theme_classic()
# 
# # Values of b to test to obtain confidence interval. For the real movement data
# # sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# # symmetric interval here.
# bvec = seq(PLB.bMLE.CI.b - 0.5, PLB.bMLE.CI.b + 0.5, 0.00001) 
# PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
# for(i in 1:length(bvec)){
#   PLB.LLvals[i] = negLL.PLB(bvec[i], x=CI.input$Area, n=length(CI.input$Area), xmin=CI.input$min.Area,
#                             xmax=CI.input$max.Area, sumlogx=CI.input$sum.log.Area)   
# }
# critVal = PLB.minNegLL.CI.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
# bIn95 = bvec[ PLB.LLvals < critVal ]
# CIbIn95 <- c(min(bIn95), max(bIn95))
# 
# # MLE Flat Rock (FR) Area ===============================================================================
# 
# 
# # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# # as a starting point for nlm for MLE of b for PLB model.
# 
# PLB.return.FR <- mle_b(Site== "Flat Rock", x=FR.input$Area, log_x=FR.input$log.Area, sum_log_x=FR.input$sum.log.Area,
#                        x_min=FR.input$min.Area, x_max=FR.input$max.Area)
# PLB.bMLE.FR.b <- PLB.return.FR[[1]] 
# PLB.minLL.FR.b <- PLB.return.FR[[2]]
# 
# # plot and find 95% confidence intervals for MLE method.
# PLB.minNegLL.FR.b <- PLB.minLL.FR.b$minimum
# x <- FR.input$Area
# FRx.PLB = seq(min(FR.input$Area), max(FR.input$Area), length=1000) # x values to plot PLB. Note
# # that these encompass the data, and are not based
# # on the binning (in MEE Figure 6 the line starts as
# # min(x), not the first bin.
# FRy.PLB = (1 - pPLB(x = FRx.PLB, b = PLB.bMLE.FR.b, xmin = min(FRx.PLB),
#                     xmax = max(FRx.PLB))) * length(FR.input$Area)
# spectra.text <- as.character(round(PLB.bMLE.FR.b, 2))
# FRb_plot <- ggplot() +
#   geom_point(aes(x = (sort(FR.input$Area, decreasing=TRUE)), y = (1:length(FR.input$Area))),
#              color = "cadetblue", size = 2, alpha = 0.3) +
#   xlab(expression(paste("Colony area, ", italic("x"), " (cm^2)"))) +
#   ylab(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    "))) +
#   scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
#                      limits = c(0.25, max(length(BR.input$Area),
#                                           length(CI.input$Area),
#                                           length(FR.input$Area),
#                                           length(FLI.input$Area),
#                                           length(IGS.input$Area),
#                                           length(JRFT.input$Area),
#                                           length(JRN.input$Area),
#                                           length(LEI.input$Area),
#                                           length(LMI.input$Area),
#                                           length(LL.input$Area),
#                                           length(Mud.input$Area),
#                                           length(NSI.input$Area),
#                                           length(NWSI.input$Area),
#                                           length(SSI.input$Area),
#                                           length(SWSI.input$Area),
#                                           length(Ten.input$Area),
#                                           length(WR.input$Area),
#                                           length(Wool.input$Area)))) +
#   scale_x_continuous(trans = 'log10', breaks = c(0,1,5,10,100,1000,10000),
#                      limits = c(min(BR.input$min.Area,
#                                     CI.input$min.Area,
#                                     FR.input$min.Area,
#                                     FLI.input$min.Area,
#                                     IGS.input$min.Area,
#                                     JRFT.input$min.Area,
#                                     JRN.input$min.Area,
#                                     LEI.input$min.Area,
#                                     LMI.input$min.Area,
#                                     LL.input$min.Area,
#                                     Mud.input$min.Area,
#                                     NSI.input$min.Area,
#                                     NWSI.input$min.Area,
#                                     SSI.input$min.Area,
#                                     SWSI.input$min.Area,
#                                     Ten.input$min.Area,
#                                     WR.input$min.Area,
#                                     Wool.input$min.Area), 
#                                 max(BR.input$max.Area,
#                                     CI.input$max.Area,
#                                     FR.input$max.Area,
#                                     FLI.input$max.Area,
#                                     IGS.input$max.Area,
#                                     JRFT.input$max.Area,
#                                     JRN.input$max.Area,
#                                     LEI.input$max.Area,
#                                     LMI.input$max.Area,
#                                     LL.input$max.Area,
#                                     Mud.input$max.Area,
#                                     NSI.input$max.Area,
#                                     NWSI.input$max.Area,
#                                     SSI.input$max.Area,
#                                     SWSI.input$max.Area,
#                                     Ten.input$max.Area,
#                                     WR.input$max.Area,
#                                     Wool.input$max.Area))) +
#   geom_line(aes(x = FRx.PLB, y = FRy.PLB), col = 'black', lwd = 1) +
#   annotate("text", x=5, y=10, label="Flat Rock") +
#   annotate("text", x=5, y=3, label = expression(paste(italic("b = "), -1.02))) +
#   theme_classic()
# 
# # Values of b to test to obtain confidence interval. For the real movement data
# # sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# # symmetric interval here.
# bvec = seq(PLB.bMLE.FR.b - 0.5, PLB.bMLE.FR.b + 0.5, 0.00001) 
# PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
# for(i in 1:length(bvec)){
#   PLB.LLvals[i] = negLL.PLB(bvec[i], x=FR.input$Area, n=length(FR.input$Area), xmin=FR.input$min.Area,
#                             xmax=FR.input$max.Area, sumlogx=FR.input$sum.log.Area)   
# }
# critVal = PLB.minNegLL.FR.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
# bIn95 = bvec[ PLB.LLvals < critVal ]
# FRbIn95 <- c(min(bIn95), max(bIn95))
# 
# # MLE Flinders (FLI) Area ===============================================================================
# # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# # as a starting point for nlm for MLE of b for PLB model.
# PLB.return.FLI <- mle_b(Site== "Flinders", x=FLI.input$Area, log_x=FLI.input$log.Area, sum_log_x=FLI.input$sum.log.Area,
#                         x_min=FLI.input$min.Area, x_max=FLI.input$max.Area)
# PLB.bMLE.FLI.b <- PLB.return.FLI[[1]] 
# PLB.minLL.FLI.b <- PLB.return.FLI[[2]]
# 
# # plot and find 95% confidence intervals for MLE method.
# PLB.minNegLL.FLI.b <- PLB.minLL.FLI.b$minimum
# x <- FLI.input$Area
# FLIx.PLB = seq(min(FLI.input$Area), max(FLI.input$Area), length=1000) # x values to plot PLB. Note
# # that these encompass the data, and are not based
# # on the binning (in MEE Figure 6 the line starts as
# # min(x), not the first bin.
# FLIy.PLB = (1 - pPLB(x = FLIx.PLB, b = PLB.bMLE.FLI.b, xmin = min(FLIx.PLB),
#                      xmax = max(FLIx.PLB))) * length(FLI.input$Area)
# spectra.text <- as.character(round(PLB.bMLE.FLI.b, 2))
# FLIb_plot <- ggplot() +
#   geom_point(aes(x = (sort(FLI.input$Area, decreasing=TRUE)), y = (1:length(FLI.input$Area))),
#              color = "darkolivegreen4", size = 2, alpha = 0.3) +
#   xlab(expression(paste("Colony area, ", italic("x"), " (cm^2)"))) +
#   ylab(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    "))) +
#   scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
#                      limits = c(0.25, max(length(BR.input$Area),
#                                           length(CI.input$Area),
#                                           length(FR.input$Area),
#                                           length(FLI.input$Area),
#                                           length(IGS.input$Area),
#                                           length(JRFT.input$Area),
#                                           length(JRN.input$Area),
#                                           length(LEI.input$Area),
#                                           length(LMI.input$Area),
#                                           length(LL.input$Area),
#                                           length(Mud.input$Area),
#                                           length(NSI.input$Area),
#                                           length(NWSI.input$Area),
#                                           length(SSI.input$Area),
#                                           length(SWSI.input$Area),
#                                           length(Ten.input$Area),
#                                           length(WR.input$Area),
#                                           length(Wool.input$Area)))) +
#   scale_x_continuous(trans = 'log10', breaks = c(0,1,5,10,100,1000,10000),
#                      limits = c(min(BR.input$min.Area,
#                                     CI.input$min.Area,
#                                     FR.input$min.Area,
#                                     FLI.input$min.Area,
#                                     IGS.input$min.Area,
#                                     JRFT.input$min.Area,
#                                     JRN.input$min.Area,
#                                     LEI.input$min.Area,
#                                     LMI.input$min.Area,
#                                     LL.input$min.Area,
#                                     Mud.input$min.Area,
#                                     NSI.input$min.Area,
#                                     NWSI.input$min.Area,
#                                     SSI.input$min.Area,
#                                     SWSI.input$min.Area,
#                                     Ten.input$min.Area,
#                                     WR.input$min.Area,
#                                     Wool.input$min.Area), 
#                                 max(BR.input$max.Area,
#                                     CI.input$max.Area,
#                                     FR.input$max.Area,
#                                     FLI.input$max.Area,
#                                     IGS.input$max.Area,
#                                     JRFT.input$max.Area,
#                                     JRN.input$max.Area,
#                                     LEI.input$max.Area,
#                                     LMI.input$max.Area,
#                                     LL.input$max.Area,
#                                     Mud.input$max.Area,
#                                     NSI.input$max.Area,
#                                     NWSI.input$max.Area,
#                                     SSI.input$max.Area,
#                                     SWSI.input$max.Area,
#                                     Ten.input$max.Area,
#                                     WR.input$max.Area,
#                                     Wool.input$max.Area))) +
#   geom_line(aes(x = FLIx.PLB, y = FLIy.PLB), col = 'black', lwd = 1) +
#   annotate("text", x=5, y=10, label="Flinders Reef") +
#   annotate("text", x=5, y=3, label = expression(paste(italic("b = "), -1.08))) +
#   theme_classic()
# 
# # Values of b to test to obtain confidence interval. For the real movement data
# # sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# # symmetric interval here.
# bvec = seq(PLB.bMLE.FLI.b - 0.5, PLB.bMLE.FLI.b + 0.5, 0.00001) 
# PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
# for(i in 1:length(bvec)){
#   PLB.LLvals[i] = negLL.PLB(bvec[i], x=FLI.input$Area, n=length(FLI.input$Area), xmin=FLI.input$min.Area,
#                             xmax=FLI.input$max.Area, sumlogx=FLI.input$sum.log.Area)   
# }
# critVal = PLB.minNegLL.FLI.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
# bIn95 = bvec[ PLB.LLvals < critVal ]
# FLIbIn95 <- c(min(bIn95), max(bIn95))
# 
# 
# # MLE Inner Gneering Shoals (IGS) Area ===============================================================================
# # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# # as a starting point for nlm for MLE of b for PLB model.
# PLB.return.IGS <- mle_b(Site== "Inner Gneering Shoals", x=IGS.input$Area, log_x=IGS.input$log.Area, sum_log_x=IGS.input$sum.log.Area,
#                         x_min=IGS.input$min.Area, x_max=IGS.input$max.Area)
# PLB.bMLE.IGS.b <- PLB.return.IGS[[1]] 
# PLB.minLL.IGS.b <- PLB.return.IGS[[2]]
# 
# # plot and find 95% confidence intervals for MLE method.
# PLB.minNegLL.IGS.b <- PLB.minLL.IGS.b$minimum
# x <- IGS.input$Area
# IGSx.PLB = seq(min(IGS.input$Area), max(IGS.input$Area), length=1000) # x values to plot PLB. Note
# # that these encompass the data, and are not based
# # on the binning (in MEE Figure 6 the line starts as
# # min(x), not the first bin.
# IGSy.PLB = (1 - pPLB(x = IGSx.PLB, b = PLB.bMLE.IGS.b, xmin = min(IGSx.PLB),
#                      xmax = max(IGSx.PLB))) * length(IGS.input$Area)
# spectra.text <- as.character(round(PLB.bMLE.IGS.b, 2))
# IGSb_plot <- ggplot() +
#   geom_point(aes(x = (sort(IGS.input$Area, decreasing=TRUE)), y = (1:length(IGS.input$Area))),
#              color = "cadetblue", size = 2, alpha = 0.3) +
#   xlab(expression(paste("Colony area, ", italic("x"), " (cm^2)"))) +
#   ylab(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    "))) +
#   scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
#                      limits = c(0.25, max(length(BR.input$Area),
#                                           length(CI.input$Area),
#                                           length(FR.input$Area),
#                                           length(FLI.input$Area),
#                                           length(IGS.input$Area),
#                                           length(JRFT.input$Area),
#                                           length(JRN.input$Area),
#                                           length(LEI.input$Area),
#                                           length(LMI.input$Area),
#                                           length(LL.input$Area),
#                                           length(Mud.input$Area),
#                                           length(NSI.input$Area),
#                                           length(NWSI.input$Area),
#                                           length(SSI.input$Area),
#                                           length(SWSI.input$Area),
#                                           length(Ten.input$Area),
#                                           length(WR.input$Area),
#                                           length(Wool.input$Area)))) +
#   scale_x_continuous(trans = 'log10', breaks = c(0,1,5,10,100,1000,10000),
#                      limits = c(min(BR.input$min.Area,
#                                     CI.input$min.Area,
#                                     FR.input$min.Area,
#                                     FLI.input$min.Area,
#                                     IGS.input$min.Area,
#                                     JRFT.input$min.Area,
#                                     JRN.input$min.Area,
#                                     LEI.input$min.Area,
#                                     LMI.input$min.Area,
#                                     LL.input$min.Area,
#                                     Mud.input$min.Area,
#                                     NSI.input$min.Area,
#                                     NWSI.input$min.Area,
#                                     SSI.input$min.Area,
#                                     SWSI.input$min.Area,
#                                     Ten.input$min.Area,
#                                     WR.input$min.Area,
#                                     Wool.input$min.Area), 
#                                 max(BR.input$max.Area,
#                                     CI.input$max.Area,
#                                     FR.input$max.Area,
#                                     FLI.input$max.Area,
#                                     IGS.input$max.Area,
#                                     JRFT.input$max.Area,
#                                     JRN.input$max.Area,
#                                     LEI.input$max.Area,
#                                     LMI.input$max.Area,
#                                     LL.input$max.Area,
#                                     Mud.input$max.Area,
#                                     NSI.input$max.Area,
#                                     NWSI.input$max.Area,
#                                     SSI.input$max.Area,
#                                     SWSI.input$max.Area,
#                                     Ten.input$max.Area,
#                                     WR.input$max.Area,
#                                     Wool.input$max.Area))) +
#   geom_line(aes(x = IGSx.PLB, y = IGSy.PLB), col = 'black', lwd = 1) +
#   annotate("text", x=10, y=10, label="Inner Gneering Shoals") +
#   annotate("text", x=5, y=3, label = expression(paste(italic("b = "), -1.04))) +
#   theme_classic()
# 
# # Values of b to test to obtain confidence interval. For the real movement data
# # sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# # symmetric interval here.
# bvec = seq(PLB.bMLE.IGS.b - 0.5, PLB.bMLE.IGS.b + 0.5, 0.00001) 
# PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
# for(i in 1:length(bvec)){
#   PLB.LLvals[i] = negLL.PLB(bvec[i], x=IGS.input$Area, n=length(IGS.input$Area), xmin=IGS.input$min.Area,
#                             xmax=IGS.input$max.Area, sumlogx=IGS.input$sum.log.Area)   
# }
# critVal = PLB.minNegLL.IGS.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
# bIn95 = bvec[ PLB.LLvals < critVal ]
# IGSbIn95 <- c(min(bIn95), max(bIn95))
# 
# # MLE Julian Rock False Trench (JRFT) Area ===============================================================================
# # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# # as a starting point for nlm for MLE of b for PLB model.
# PLB.return.JRFT <- mle_b(Site== "Julian Rock False Trench", x=JRFT.input$Area, log_x=JRFT.input$log.Area, sum_log_x=JRFT.input$sum.log.Area,
#                          x_min=JRFT.input$min.Area, x_max=JRFT.input$max.Area)
# PLB.bMLE.JRFT.b <- PLB.return.JRFT[[1]] 
# PLB.minLL.JRFT.b <- PLB.return.JRFT[[2]]
# 
# # plot and find 95% confidence intervals for MLE method.
# PLB.minNegLL.JRFT.b <- PLB.minLL.JRFT.b$minimum
# x <- JRFT.input$Area
# JRFTx.PLB = seq(min(JRFT.input$Area), max(JRFT.input$Area), length=1000) # x values to plot PLB. Note
# # that these encompass the data, and are not based
# # on the binning (in MEE Figure 6 the line starts as
# # min(x), not the first bin.
# JRFTy.PLB = (1 - pPLB(x = JRFTx.PLB, b = PLB.bMLE.JRFT.b, xmin = min(JRFTx.PLB),
#                       xmax = max(JRFTx.PLB))) * length(JRFT.input$Area)
# spectra.text <- as.character(round(PLB.bMLE.JRFT.b, 2))
# JRFTb_plot <- ggplot() +
#   geom_point(aes(x = (sort(JRFT.input$Area, decreasing=TRUE)), y = (1:length(JRFT.input$Area))),
#              color = "cadetblue", size = 2, alpha = 0.3) +
#   xlab(expression(paste("Colony area, ", italic("x"), " (cm^2)"))) +
#   ylab(expression(paste("Number of colonies", " ">=" ", italic("x"), "    "))) +
#   scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
#                      limits = c(0.25, max(length(BR.input$Area),
#                                           length(CI.input$Area),
#                                           length(FR.input$Area),
#                                           length(FLI.input$Area),
#                                           length(IGS.input$Area),
#                                           length(JRFT.input$Area),
#                                           length(JRN.input$Area),
#                                           length(LEI.input$Area),
#                                           length(LMI.input$Area),
#                                           length(LL.input$Area),
#                                           length(Mud.input$Area),
#                                           length(NSI.input$Area),
#                                           length(NWSI.input$Area),
#                                           length(SSI.input$Area),
#                                           length(SWSI.input$Area),
#                                           length(Ten.input$Area),
#                                           length(WR.input$Area),
#                                           length(Wool.input$Area)))) +
#   scale_x_continuous(trans = 'log10', breaks = c(0,1,5,10,100,1000,10000),
#                      limits = c(min(BR.input$min.Area,
#                                     CI.input$min.Area,
#                                     FR.input$min.Area,
#                                     FLI.input$min.Area,
#                                     IGS.input$min.Area,
#                                     JRFT.input$min.Area,
#                                     JRN.input$min.Area,
#                                     LEI.input$min.Area,
#                                     LMI.input$min.Area,
#                                     LL.input$min.Area,
#                                     Mud.input$min.Area,
#                                     NSI.input$min.Area,
#                                     NWSI.input$min.Area,
#                                     SSI.input$min.Area,
#                                     SWSI.input$min.Area,
#                                     Ten.input$min.Area,
#                                     WR.input$min.Area,
#                                     Wool.input$min.Area), 
#                                 max(BR.input$max.Area,
#                                     CI.input$max.Area,
#                                     FR.input$max.Area,
#                                     FLI.input$max.Area,
#                                     IGS.input$max.Area,
#                                     JRFT.input$max.Area,
#                                     JRN.input$max.Area,
#                                     LEI.input$max.Area,
#                                     LMI.input$max.Area,
#                                     LL.input$max.Area,
#                                     Mud.input$max.Area,
#                                     NSI.input$max.Area,
#                                     NWSI.input$max.Area,
#                                     SSI.input$max.Area,
#                                     SWSI.input$max.Area,
#                                     Ten.input$max.Area,
#                                     WR.input$max.Area,
#                                     Wool.input$max.Area))) +
#   geom_line(aes(x = JRFTx.PLB, y = JRFTy.PLB), col = 'black', lwd = 1) +
#   annotate("text", x=10, y=10, label="Julian Rock False Trench") +
#   annotate("text", x=5, y=3, label = expression(paste(italic("b = "), -0.953))) +
#   theme_classic()
# 
# # Values of b to test to obtain confidence interval. For the real movement data
# # sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# # symmetric interval here.
# bvec = seq(PLB.bMLE.JRFT.b - 0.5, PLB.bMLE.JRFT.b + 0.5, 0.00001) 
# PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
# for(i in 1:length(bvec)){
#   PLB.LLvals[i] = negLL.PLB(bvec[i], x=JRFT.input$Area, n=length(JRFT.input$Area), xmin=JRFT.input$min.Area,
#                             xmax=JRFT.input$max.Area, sumlogx=JRFT.input$sum.log.Area)   
# }
# critVal = PLB.minNegLL.JRFT.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
# bIn95 = bvec[ PLB.LLvals < critVal ]
# JRFTbIn95 <- c(min(bIn95), max(bIn95))
# 
# # MLE Julian Rock Nursery (JRN) Area ===============================================================================
# # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# # as a starting point for nlm for MLE of b for PLB model.
# PLB.return.JRN <- mle_b(Site== "Julian Rock Nursery", x=JRN.input$Area, log_x=JRN.input$log.Area, sum_log_x=JRN.input$sum.log.Area,
#                         x_min=JRN.input$min.Area, x_max=JRN.input$max.Area)
# PLB.bMLE.JRN.b <- PLB.return.JRN[[1]] 
# PLB.minLL.JRN.b <- PLB.return.JRN[[2]]
# 
# # plot and find 95% confidence intervals for MLE method.
# PLB.minNegLL.JRN.b <- PLB.minLL.JRN.b$minimum
# x <- JRN.input$Area
# JRNx.PLB = seq(min(JRN.input$Area), max(JRN.input$Area), length=1000) # x values to plot PLB. Note
# # that these encompass the data, and are not based
# # on the binning (in MEE Figure 6 the line starts as
# # min(x), not the first bin.
# JRNy.PLB = (1 - pPLB(x = JRNx.PLB, b = PLB.bMLE.JRN.b, xmin = min(JRNx.PLB),
#                      xmax = max(JRNx.PLB))) * length(JRN.input$Area)
# spectra.text <- as.character(round(PLB.bMLE.JRN.b, 2))
# JRNb_plot <- ggplot() +
#   geom_point(aes(x = (sort(JRN.input$Area, decreasing=TRUE)), y = (1:length(JRN.input$Area))),
#              color = "cadetblue", size = 2, alpha = 0.3) +
#   xlab(expression(paste("Colony area, ", italic("x"), " (cm^2)"))) +
#   ylab(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    "))) +
#   scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
#                      limits = c(0.25, max(length(BR.input$Area),
#                                           length(CI.input$Area),
#                                           length(FR.input$Area),
#                                           length(FLI.input$Area),
#                                           length(IGS.input$Area),
#                                           length(JRFT.input$Area),
#                                           length(JRN.input$Area),
#                                           length(LEI.input$Area),
#                                           length(LMI.input$Area),
#                                           length(LL.input$Area),
#                                           length(Mud.input$Area),
#                                           length(NSI.input$Area),
#                                           length(NWSI.input$Area),
#                                           length(SSI.input$Area),
#                                           length(SWSI.input$Area),
#                                           length(Ten.input$Area),
#                                           length(WR.input$Area),
#                                           length(Wool.input$Area)))) +
#   scale_x_continuous(trans = 'log10', breaks = c(0,1,5,10,100,1000,10000),
#                      limits = c(min(BR.input$min.Area,
#                                     CI.input$min.Area,
#                                     FR.input$min.Area,
#                                     FLI.input$min.Area,
#                                     IGS.input$min.Area,
#                                     JRFT.input$min.Area,
#                                     JRN.input$min.Area,
#                                     LEI.input$min.Area,
#                                     LMI.input$min.Area,
#                                     LL.input$min.Area,
#                                     Mud.input$min.Area,
#                                     NSI.input$min.Area,
#                                     NWSI.input$min.Area,
#                                     SSI.input$min.Area,
#                                     SWSI.input$min.Area,
#                                     Ten.input$min.Area,
#                                     WR.input$min.Area,
#                                     Wool.input$min.Area), 
#                                 max(BR.input$max.Area,
#                                     CI.input$max.Area,
#                                     FR.input$max.Area,
#                                     FLI.input$max.Area,
#                                     IGS.input$max.Area,
#                                     JRFT.input$max.Area,
#                                     JRN.input$max.Area,
#                                     LEI.input$max.Area,
#                                     LMI.input$max.Area,
#                                     LL.input$max.Area,
#                                     Mud.input$max.Area,
#                                     NSI.input$max.Area,
#                                     NWSI.input$max.Area,
#                                     SSI.input$max.Area,
#                                     SWSI.input$max.Area,
#                                     Ten.input$max.Area,
#                                     WR.input$max.Area,
#                                     Wool.input$max.Area))) +
#   geom_line(aes(x = JRNx.PLB, y = JRNy.PLB), col = 'black', lwd = 1) +
#   annotate("text", x=5, y=10, label="Julian Rock Nursery") +
#   annotate("text", x=5, y=3, label = expression(paste(italic("b = "), -1.05))) +
#   theme_classic()
# 
# # Values of b to test to obtain confidence interval. For the real movement data
# # sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# # symmetric interval here.
# bvec = seq(PLB.bMLE.JRN.b - 0.5, PLB.bMLE.JRN.b + 0.5, 0.00001) 
# PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
# for(i in 1:length(bvec)){
#   PLB.LLvals[i] = negLL.PLB(bvec[i], x=JRN.input$Area, n=length(JRN.input$Area), xmin=JRN.input$min.Area,
#                             xmax=JRN.input$max.Area, sumlogx=JRN.input$sum.log.Area)   
# }
# critVal = PLB.minNegLL.JRN.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
# bIn95 = bvec[ PLB.LLvals < critVal ]
# JRNbIn95 <- c(min(bIn95), max(bIn95))
# 
# # MLE Lady Elliot Island (LEI) Area ===============================================================================
# 
# # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# # as a starting point for nlm for MLE of b for PLB model.
# PLB.return.LEI <- mle_b(Site== "Lady Elliot Island", x=LEI.input$Area, log_x=LEI.input$log.Area, sum_log_x=LEI.input$sum.log.Area,
#                         x_min=LEI.input$min.Area, x_max=LEI.input$max.Area)
# PLB.bMLE.LEI.b <- PLB.return.LEI[[1]] 
# PLB.minLL.LEI.b <- PLB.return.LEI[[2]]
# 
# # plot and find 95% confidence intervals for MLE method.
# PLB.minNegLL.LEI.b <- PLB.minLL.LEI.b$minimum
# x <- LEI.input$Area
# LEIx.PLB = seq(min(LEI.input$Area), max(LEI.input$Area), length=1000) # x values to plot PLB. Note
# # that these encompass the data, and are not based
# # on the binning (in MEE Figure 6 the line starts as
# # min(x), not the first bin.
# LEIy.PLB = (1 - pPLB(x = LEIx.PLB, b = PLB.bMLE.LEI.b, xmin = min(LEIx.PLB),
#                      xmax = max(LEIx.PLB))) * length(LEI.input$Area)
# spectra.text <- as.character(round(PLB.bMLE.LEI.b, 2))
# LEIb_plot <- ggplot() +
#   geom_point(aes(x = (sort(LEI.input$Area, decreasing=TRUE)), y = (1:length(LEI.input$Area))), 
#              color = "deeppink3", size = 2, alpha = 0.3) +
#   xlab(expression(paste("Colony area, ", italic("x"), " (cm^2)"))) +
#   ylab(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    ")))+
#   scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
#                      limits = c(0.25, max(length(BR.input$Area),
#                                           length(CI.input$Area),
#                                           length(FR.input$Area),
#                                           length(FLI.input$Area),
#                                           length(IGS.input$Area),
#                                           length(JRFT.input$Area),
#                                           length(JRN.input$Area),
#                                           length(LEI.input$Area),
#                                           length(LMI.input$Area),
#                                           length(LL.input$Area),
#                                           length(Mud.input$Area),
#                                           length(NSI.input$Area),
#                                           length(NWSI.input$Area),
#                                           length(SSI.input$Area),
#                                           length(SWSI.input$Area),
#                                           length(Ten.input$Area),
#                                           length(WR.input$Area),
#                                           length(Wool.input$Area)))) +
#   scale_x_continuous(trans = 'log10', breaks = c(1,5,10,100,1000,10000),
#                      limits = c(min(BR.input$min.Area,
#                                     CI.input$min.Area,
#                                     FR.input$min.Area,
#                                     FLI.input$min.Area,
#                                     IGS.input$min.Area,
#                                     JRFT.input$min.Area,
#                                     JRN.input$min.Area,
#                                     LEI.input$min.Area,
#                                     LMI.input$min.Area,
#                                     LL.input$min.Area,
#                                     Mud.input$min.Area,
#                                     NSI.input$min.Area,
#                                     NWSI.input$min.Area,
#                                     SSI.input$min.Area,
#                                     SWSI.input$min.Area,
#                                     Ten.input$min.Area,
#                                     WR.input$min.Area,
#                                     Wool.input$min.Area), 
#                                 max(BR.input$max.Area,
#                                     CI.input$max.Area,
#                                     FR.input$max.Area,
#                                     FLI.input$max.Area,
#                                     IGS.input$max.Area,
#                                     JRFT.input$max.Area,
#                                     JRN.input$max.Area,
#                                     LEI.input$max.Area,
#                                     LMI.input$max.Area,
#                                     LL.input$max.Area,
#                                     Mud.input$max.Area,
#                                     NSI.input$max.Area,
#                                     NWSI.input$max.Area,
#                                     SSI.input$max.Area,
#                                     SWSI.input$max.Area,
#                                     Ten.input$max.Area,
#                                     WR.input$max.Area,
#                                     Wool.input$max.Area))) +
#   geom_line(aes(x = LEIx.PLB, y = LEIy.PLB), col = 'black', lwd = 1) +
#   annotate("text", x=5, y=10, label="Lady Elliot Island") +
#   annotate("text", x=5, y=3, label = expression(paste(italic("b = "), -0.996))) +
#   theme_classic()
# 
# # Values of b to test to obtain confidence interval. For the real movement data
# # sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# # symmetric interval here.
# bvec = seq(PLB.bMLE.LEI.b - 0.5, PLB.bMLE.LEI.b + 0.5, 0.00001) 
# PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
# for(i in 1:length(bvec)){
#   PLB.LLvals[i] = negLL.PLB(bvec[i], x=LL.input$Area, n=length(LL.input$Area), xmin=LL.input$min.Area,
#                             xmax=LL.input$max.Area, sumlogx=LL.input$sum.log.Area)   
# }
# critVal = PLB.minNegLL.LEI.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
# bIn95 = bvec[ PLB.LLvals < critVal ]
# lobIn95 <- c(min(bIn95), max(bIn95))
# 
# # MLE Lady Musgrave (LMI) Area ===============================================================================
# 
# # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# # as a starting point for nlm for MLE of b for PLB model.
# PLB.return.LMI <- mle_b(Site== "Lady Musgrave", x=LMI.input$Area, log_x=LMI.input$log.Area, sum_log_x=LMI.input$sum.log.Area,
#                        x_min=LMI.input$min.Area, x_max=LMI.input$max.Area)
# PLB.bMLE.LMI.b <- PLB.return.LMI[[1]] 
# PLB.minLL.LMI.b <- PLB.return.LMI[[2]]
# 
# # plot and find 95% confidence intervals for MLE method.
# PLB.minNegLL.LMI.b <- PLB.minLL.LMI.b$minimum
# x <- LMI.input$Area
# LMIx.PLB = seq(min(LMI.input$Area), max(LMI.input$Area), length=1000) # x values to plot PLB. Note
# # that these encompass the data, and are not based
# # on the binning (in MEE Figure 6 the line starts as
# # min(x), not the first bin.
# LMIy.PLB = (1 - pPLB(x = LMIx.PLB, b = PLB.bMLE.LMI.b, xmin = min(LMIx.PLB),
#                     xmax = max(LMIx.PLB))) * length(LMI.input$Area)
# spectra.text <- as.character(round(PLB.bMLE.LMI.b, 2))
# LMIb_plot <- ggplot() +
#   geom_point(aes(x = (sort(LMI.input$Area, decreasing=TRUE)), y = (1:length(LMI.input$Area))), 
#              color = "deeppink3", size = 2, alpha = 0.3) +
#   xlab(expression(paste("Colony area, ", italic("x"), " (cm^2)"))) +
#   ylab(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    "))) +
#   scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
#                      limits = c(0.25, max(length(BR.input$Area),
#                                           length(CI.input$Area),
#                                           length(FR.input$Area),
#                                           length(FLI.input$Area),
#                                           length(IGS.input$Area),
#                                           length(JRFT.input$Area),
#                                           length(JRN.input$Area),
#                                           length(LEI.input$Area),
#                                           length(LMI.input$Area),
#                                           length(LL.input$Area),
#                                           length(Mud.input$Area),
#                                           length(NSI.input$Area),
#                                           length(NWSI.input$Area),
#                                           length(SSI.input$Area),
#                                           length(SWSI.input$Area),
#                                           length(Ten.input$Area),
#                                           length(WR.input$Area),
#                                           length(Wool.input$Area)))) +
#   scale_x_continuous(trans = 'log10', breaks = c(1,5,10,100,1000,10000),
#                      limits = c(min(BR.input$min.Area,
#                                     CI.input$min.Area,
#                                     FR.input$min.Area,
#                                     FLI.input$min.Area,
#                                     IGS.input$min.Area,
#                                     JRFT.input$min.Area,
#                                     JRN.input$min.Area,
#                                     LEI.input$min.Area,
#                                     LMI.input$min.Area,
#                                     LL.input$min.Area,
#                                     Mud.input$min.Area,
#                                     NSI.input$min.Area,
#                                     NWSI.input$min.Area,
#                                     SSI.input$min.Area,
#                                     SWSI.input$min.Area,
#                                     Ten.input$min.Area,
#                                     WR.input$min.Area,
#                                     Wool.input$min.Area), 
#                                 max(BR.input$max.Area,
#                                     CI.input$max.Area,
#                                     FR.input$max.Area,
#                                     FLI.input$max.Area,
#                                     IGS.input$max.Area,
#                                     JRFT.input$max.Area,
#                                     JRN.input$max.Area,
#                                     LEI.input$max.Area,
#                                     LMI.input$max.Area,
#                                     LL.input$max.Area,
#                                     Mud.input$max.Area,
#                                     NSI.input$max.Area,
#                                     NWSI.input$max.Area,
#                                     SSI.input$max.Area,
#                                     SWSI.input$max.Area,
#                                     Ten.input$max.Area,
#                                     WR.input$max.Area,
#                                     Wool.input$max.Area))) +
#   geom_line(aes(x = LMIx.PLB, y = LMIy.PLB), col = 'black', lwd = 1) +
#   annotate("text", x=10, y=10, label="Lady Musgrave Island") +
#   annotate("text", x=5, y=3, label = expression(paste(italic("  b = "), -1.10))) +
#   theme_classic()
# 
# # Values of b to test to obtain confidence interval. For the real movement data
# # sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# # symmetric interval here.
# bvec = seq(PLB.bMLE.LMI.b - 0.5, PLB.bMLE.LMI.b + 0.5, 0.00001) 
# PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
# for(i in 1:length(bvec)){
#   PLB.LLvals[i] = negLL.PLB(bvec[i], x=LMI.input$Area, n=length(LMI.input$Area), xmin=LMI.input$min.Area,
#                             xmax=LMI.input$max.Area, sumlogx=LMI.input$sum.log.Area)   
# }
# critVal = PLB.minNegLL.LMI.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
# bIn95 = bvec[ PLB.LLvals < critVal ]
# LMIbIn95 <- c(min(bIn95), max(bIn95))
# 
# # MLE Libbys Lair (LL) Area ===============================================================================
# # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# # as a starting point for nlm for MLE of b for PLB model.
# PLB.return.LL <- mle_b(Site== "Libbys Lair", x=LL.input$Area, log_x=LL.input$log.Area, sum_log_x=LL.input$sum.log.Area,
#                        x_min=LL.input$min.Area, x_max=LL.input$max.Area)
# PLB.bMLE.LL.b <- PLB.return.LL[[1]] 
# PLB.minLL.LL.b <- PLB.return.LL[[2]]
# 
# # plot and find 95% confidence intervals for MLE method.
# PLB.minNegLL.LL.b <- PLB.minLL.LL.b$minimum
# x <- LL.input$Area
# LLx.PLB = seq(min(LL.input$Area), max(LL.input$Area), length=1000) # x values to plot PLB. Note
# # that these encompass the data, and are not based
# # on the binning (in MEE Figure 6 the line starts as
# # min(x), not the first bin.
# LLy.PLB = (1 - pPLB(x = LLx.PLB, b = PLB.bMLE.LL.b, xmin = min(LLx.PLB),
#                     xmax = max(LLx.PLB))) * length(LL.input$Area)
# spectra.text <- as.character(round(PLB.bMLE.LL.b, 2))
# LLb_plot <- ggplot() +
#   geom_point(aes(x = (sort(LL.input$Area, decreasing=TRUE)), y = (1:length(LL.input$Area))),
#              color = "deeppink3", size = 2, alpha = 0.3) +
#   xlab(expression(paste("Colony area, ", italic("x"), " (cm^2)"))) +
#   ylab(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    "))) +
#   scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
#                      limits = c(0.25, max(length(BR.input$Area),
#                                           length(CI.input$Area),
#                                           length(FR.input$Area),
#                                           length(FLI.input$Area),
#                                           length(IGS.input$Area),
#                                           length(JRFT.input$Area),
#                                           length(JRN.input$Area),
#                                           length(LEI.input$Area),
#                                           length(LMI.input$Area),
#                                           length(LL.input$Area),
#                                           length(Mud.input$Area),
#                                           length(NSI.input$Area),
#                                           length(NWSI.input$Area),
#                                           length(SSI.input$Area),
#                                           length(SWSI.input$Area),
#                                           length(Ten.input$Area),
#                                           length(WR.input$Area),
#                                           length(Wool.input$Area)))) +
#   scale_x_continuous(trans = 'log10', breaks = c(0,1,5,10,100,1000,10000),
#                      limits = c(min(BR.input$min.Area,
#                                     CI.input$min.Area,
#                                     FR.input$min.Area,
#                                     FLI.input$min.Area,
#                                     IGS.input$min.Area,
#                                     JRFT.input$min.Area,
#                                     JRN.input$min.Area,
#                                     LEI.input$min.Area,
#                                     LMI.input$min.Area,
#                                     LL.input$min.Area,
#                                     Mud.input$min.Area,
#                                     NSI.input$min.Area,
#                                     NWSI.input$min.Area,
#                                     SSI.input$min.Area,
#                                     SWSI.input$min.Area,
#                                     Ten.input$min.Area,
#                                     WR.input$min.Area,
#                                     Wool.input$min.Area), 
#                                 max(BR.input$max.Area,
#                                     CI.input$max.Area,
#                                     FR.input$max.Area,
#                                     FLI.input$max.Area,
#                                     IGS.input$max.Area,
#                                     JRFT.input$max.Area,
#                                     JRN.input$max.Area,
#                                     LEI.input$max.Area,
#                                     LMI.input$max.Area,
#                                     LL.input$max.Area,
#                                     Mud.input$max.Area,
#                                     NSI.input$max.Area,
#                                     NWSI.input$max.Area,
#                                     SSI.input$max.Area,
#                                     SWSI.input$max.Area,
#                                     Ten.input$max.Area,
#                                     WR.input$max.Area,
#                                     Wool.input$max.Area))) +
#   geom_line(aes(x = LLx.PLB, y = LLy.PLB), col = 'black', lwd = 1) +
#   annotate("text", x=5, y=10, label="Libby's Lair") +
#   annotate("text", x=5, y=3, label = expression(paste(italic("b = "), -1.09))) +
#   theme_classic()
# 
# # Values of b to test to obtain confidence interval. For the real movement data
# # sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# # symmetric interval here.
# bvec = seq(PLB.bMLE.LL.b - 0.5, PLB.bMLE.LL.b + 0.5, 0.00001) 
# PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
# for(i in 1:length(bvec)){
#   PLB.LLvals[i] = negLL.PLB(bvec[i], x=LL.input$Area, n=length(LL.input$Area), xmin=LL.input$min.Area,
#                             xmax=LL.input$max.Area, sumlogx=LL.input$sum.log.Area)   
# }
# critVal = PLB.minNegLL.LL.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
# bIn95 = bvec[ PLB.LLvals < critVal ]
# LLbIn95 <- c(min(bIn95), max(bIn95))
# 
# # MLE Mudjimba (Mud) Area ===============================================================================
# 
# # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# # as a starting point for nlm for MLE of b for PLB model.
# PLB.return.Mud <- mle_b(Site== "Mudjimba", x=Mud.input$Area, log_x=Mud.input$log.Area, sum_log_x=Mud.input$sum.log.Area,
#                         x_min=Mud.input$min.Area, x_max=Mud.input$max.Area)
# PLB.bMLE.Mud.b <- PLB.return.Mud[[1]] 
# PLB.minLL.Mud.b <- PLB.return.Mud[[2]]
# 
# # plot and find 95% confidence intervals for MLE method.
# PLB.minNegLL.Mud.b <- PLB.minLL.Mud.b$minimum
# x <- Mud.input$Area
# Mudx.PLB = seq(min(Mud.input$Area), max(Mud.input$Area), length=1000) # x values to plot PLB. Note
# # that these encompass the data, and are not based
# # on the binning (in MEE Figure 6 the line starts as
# # min(x), not the first bin.
# Mudy.PLB = (1 - pPLB(x = Mudx.PLB, b = PLB.bMLE.Mud.b, xmin = min(Mudx.PLB),
#                      xmax = max(Mudx.PLB))) * length(Mud.input$Area)
# spectra.text <- as.character(round(PLB.bMLE.Mud.b, 2))
# Mudb_plot <- ggplot() +
#   geom_point(aes(x = (sort(Mud.input$Area, decreasing=TRUE)), y = (1:length(Mud.input$Area))), 
#              color = "cadetblue", size = 2, alpha = 0.3) +
#   xlab(expression(paste("Colony area, ", italic("x"), " (cm^2)"))) +
#   ylab(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    ")))+
#   scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
#                      limits = c(0.25, max(length(BR.input$Area),
#                                           length(CI.input$Area),
#                                           length(FR.input$Area),
#                                           length(FLI.input$Area),
#                                           length(IGS.input$Area),
#                                           length(JRFT.input$Area),
#                                           length(JRN.input$Area),
#                                           length(LEI.input$Area),
#                                           length(LMI.input$Area),
#                                           length(LL.input$Area),
#                                           length(Mud.input$Area),
#                                           length(NSI.input$Area),
#                                           length(NWSI.input$Area),
#                                           length(SSI.input$Area),
#                                           length(SWSI.input$Area),
#                                           length(Ten.input$Area),
#                                           length(WR.input$Area),
#                                           length(Wool.input$Area)))) +
#   scale_x_continuous(trans = 'log10', breaks = c(1,5,10,100,1000,10000),
#                      limits = c(min(BR.input$min.Area,
#                                     CI.input$min.Area,
#                                     FR.input$min.Area,
#                                     FLI.input$min.Area,
#                                     IGS.input$min.Area,
#                                     JRFT.input$min.Area,
#                                     JRN.input$min.Area,
#                                     LEI.input$min.Area,
#                                     LMI.input$min.Area,
#                                     LL.input$min.Area,
#                                     Mud.input$min.Area,
#                                     NSI.input$min.Area,
#                                     NWSI.input$min.Area,
#                                     SSI.input$min.Area,
#                                     SWSI.input$min.Area,
#                                     Ten.input$min.Area,
#                                     WR.input$min.Area,
#                                     Wool.input$min.Area), 
#                                 max(BR.input$max.Area,
#                                     CI.input$max.Area,
#                                     FR.input$max.Area,
#                                     FLI.input$max.Area,
#                                     IGS.input$max.Area,
#                                     JRFT.input$max.Area,
#                                     JRN.input$max.Area,
#                                     LEI.input$max.Area,
#                                     LMI.input$max.Area,
#                                     LL.input$max.Area,
#                                     Mud.input$max.Area,
#                                     NSI.input$max.Area,
#                                     NWSI.input$max.Area,
#                                     SSI.input$max.Area,
#                                     SWSI.input$max.Area,
#                                     Ten.input$max.Area,
#                                     WR.input$max.Area,
#                                     Wool.input$max.Area))) +
#   geom_line(aes(x = Mudx.PLB, y = Mudy.PLB), col = 'black', lwd = 1) +
#   annotate("text", x=5, y=10, label="Mudjimba") +
#   annotate("text", x=5, y=3, label = expression(paste(italic("b = "), -1.14))) +
#   theme_classic()
# 
# # Values of b to test to obtain confidence interval. For the real movement data
# # sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# # symmetric interval here.
# bvec = seq(PLB.bMLE.Mud.b - 0.5, PLB.bMLE.Mud.b + 0.5, 0.00001) 
# PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
# for(i in 1:length(bvec)){
#   PLB.LLvals[i] = negLL.PLB(bvec[i], x=Mud.input$Area, n=length(Mud.input$Area), xmin=Mud.input$min.Area,
#                             xmax=Mud.input$max.Area, sumlogx=Mud.input$sum.log.Area)   
# }
# critVal = PLB.minNegLL.Mud.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
# bIn95 = bvec[ PLB.LLvals < critVal ]
# MudbIn95 <- c(min(bIn95), max(bIn95))
# 
# # MLE North Solitary Island (NSI) Area ===============================================================================
# 
# # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# # as a starting point for nlm for MLE of b for PLB model.
# PLB.return.NSI <- mle_b(Site== "North Solitary Island", x=NSI.input$Area, log_x=NSI.input$log.Area, sum_log_x=NSI.input$sum.log.Area,
#                         x_min=NSI.input$min.Area, x_max=NSI.input$max.Area)
# PLB.bMLE.NSI.b <- PLB.return.NSI[[1]] 
# PLB.minLL.NSI.b <- PLB.return.NSI[[2]]
# 
# # plot and find 95% confidence intervals for MLE method.
# PLB.minNegLL.NSI.b <- PLB.minLL.NSI.b$minimum
# x <- NSI.input$Area
# NSIx.PLB = seq(min(NSI.input$Area), max(NSI.input$Area), length=1000) # x values to plot PLB. Note
# # that these encompass the data, and are not based
# # on the binning (in MEE Figure 6 the line starts as
# # min(x), not the first bin.
# NSIy.PLB = (1 - pPLB(x = NSIx.PLB, b = PLB.bMLE.NSI.b, xmin = min(NSIx.PLB),
#                      xmax = max(NSIx.PLB))) * length(NSI.input$Area)
# spectra.text <- as.character(round(PLB.bMLE.NSI.b, 2))
# NSIb_plot <- ggplot() +
#   geom_point(aes(x = (sort(NSI.input$Area, decreasing=TRUE)), y = (1:length(NSI.input$Area))), 
#              color = "cadetblue", size = 2, alpha = 0.3) +
#   xlab(expression(paste("Colony area, ", italic("x"), " (cm^2)"))) +
#   ylab(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    ")))+
#   scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
#                      limits = c(0.25, max(length(BR.input$Area),
#                                           length(CI.input$Area),
#                                           length(FR.input$Area),
#                                           length(FLI.input$Area),
#                                           length(IGS.input$Area),
#                                           length(JRFT.input$Area),
#                                           length(JRN.input$Area),
#                                           length(LEI.input$Area),
#                                           length(LMI.input$Area),
#                                           length(LL.input$Area),
#                                           length(Mud.input$Area),
#                                           length(NSI.input$Area),
#                                           length(NWSI.input$Area),
#                                           length(SSI.input$Area),
#                                           length(SWSI.input$Area),
#                                           length(Ten.input$Area),
#                                           length(WR.input$Area),
#                                           length(Wool.input$Area)))) +
#   scale_x_continuous(trans = 'log10', breaks = c(1,5,10,100,1000,10000),
#                      limits = c(min(BR.input$min.Area,
#                                     CI.input$min.Area,
#                                     FR.input$min.Area,
#                                     FLI.input$min.Area,
#                                     IGS.input$min.Area,
#                                     JRFT.input$min.Area,
#                                     JRN.input$min.Area,
#                                     LEI.input$min.Area,
#                                     LMI.input$min.Area,
#                                     LL.input$min.Area,
#                                     Mud.input$min.Area,
#                                     NSI.input$min.Area,
#                                     NWSI.input$min.Area,
#                                     SSI.input$min.Area,
#                                     SWSI.input$min.Area,
#                                     Ten.input$min.Area,
#                                     WR.input$min.Area,
#                                     Wool.input$min.Area), 
#                                 max(BR.input$max.Area,
#                                     CI.input$max.Area,
#                                     FR.input$max.Area,
#                                     FLI.input$max.Area,
#                                     IGS.input$max.Area,
#                                     JRFT.input$max.Area,
#                                     JRN.input$max.Area,
#                                     LEI.input$max.Area,
#                                     LMI.input$max.Area,
#                                     LL.input$max.Area,
#                                     Mud.input$max.Area,
#                                     NSI.input$max.Area,
#                                     NWSI.input$max.Area,
#                                     SSI.input$max.Area,
#                                     SWSI.input$max.Area,
#                                     Ten.input$max.Area,
#                                     WR.input$max.Area,
#                                     Wool.input$max.Area))) +
#   geom_line(aes(x = NSIx.PLB, y = NSIy.PLB), col = 'black', lwd = 1) +
#   annotate("text", x=5, y=10, label="North Solitary Island") +
#   annotate("text", x=5, y=3, label = expression(paste(italic("b = "), -0.969))) +
#   theme_classic()
# 
# # Values of b to test to obtain confidence interval. For the real movement data
# # sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# # symmetric interval here.
# bvec = seq(PLB.bMLE.NSI.b - 0.5, PLB.bMLE.NSI.b + 0.5, 0.00001) 
# PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
# for(i in 1:length(bvec)){
#   PLB.LLvals[i] = negLL.PLB(bvec[i], x=NSI.input$Area, n=length(NSI.input$Area), xmin=NSI.input$min.Area,
#                             xmax=NSI.input$max.Area, sumlogx=NSI.input$sum.log.Area)   
# }
# critVal = PLB.minNegLL.NSI.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
# bIn95 = bvec[ PLB.LLvals < critVal ]
# NSIbIn95 <- c(min(bIn95), max(bIn95))
# 
# 
# # MLE North West Solitary Island (NWSI) Area ===============================================================================
# # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# # as a starting point for nlm for MLE of b for PLB model.
# PLB.return.NWSI <- mle_b(Site== "North West Solitary Island", x=NWSI.input$Area, log_x=NWSI.input$log.Area, sum_log_x=NWSI.input$sum.log.Area,
#                          x_min=NWSI.input$min.Area, x_max=NWSI.input$max.Area)
# PLB.bMLE.NWSI.b <- PLB.return.NWSI[[1]] 
# PLB.minLL.NWSI.b <- PLB.return.NWSI[[2]]
# 
# # plot and find 95% confidence intervals for MLE method.
# PLB.minNegLL.NWSI.b <- PLB.minLL.NWSI.b$minimum
# x <- NWSI.input$Area
# NWSIx.PLB = seq(min(NWSI.input$Area), max(NWSI.input$Area), length=1000) # x values to plot PLB. Note
# # that these encompass the data, and are not based
# # on the binning (in MEE Figure 6 the line starts as
# # min(x), not the first bin.
# NWSIy.PLB = (1 - pPLB(x = NWSIx.PLB, b = PLB.bMLE.NWSI.b, xmin = min(NWSIx.PLB),
#                       xmax = max(NWSIx.PLB))) * length(NWSI.input$Area)
# spectra.text <- as.character(round(PLB.bMLE.NWSI.b, 2))
# NWSIb_plot <- ggplot() +
#   geom_point(aes(x = (sort(NWSI.input$Area, decreasing=TRUE)), y = (1:length(NWSI.input$Area))),
#              color = "cadetblue", size = 2, alpha = 0.3) +
#   xlab(expression(paste("Colony area, ", italic("x"), " (cm^2)"))) +
#   ylab(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    "))) +
#   scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
#                      limits = c(0.25, max(length(BR.input$Area),
#                                           length(CI.input$Area),
#                                           length(FR.input$Area),
#                                           length(FLI.input$Area),
#                                           length(IGS.input$Area),
#                                           length(JRFT.input$Area),
#                                           length(JRN.input$Area),
#                                           length(LEI.input$Area),
#                                           length(LMI.input$Area),
#                                           length(LL.input$Area),
#                                           length(Mud.input$Area),
#                                           length(NSI.input$Area),
#                                           length(NWSI.input$Area),
#                                           length(SSI.input$Area),
#                                           length(SWSI.input$Area),
#                                           length(Ten.input$Area),
#                                           length(WR.input$Area),
#                                           length(Wool.input$Area)))) +
#   scale_x_continuous(trans = 'log10', breaks = c(0,1,5,10,100,1000,10000),
#                      limits = c(min(BR.input$min.Area,
#                                     CI.input$min.Area,
#                                     FR.input$min.Area,
#                                     FLI.input$min.Area,
#                                     IGS.input$min.Area,
#                                     JRFT.input$min.Area,
#                                     JRN.input$min.Area,
#                                     LEI.input$min.Area,
#                                     LMI.input$min.Area,
#                                     LL.input$min.Area,
#                                     Mud.input$min.Area,
#                                     NSI.input$min.Area,
#                                     NWSI.input$min.Area,
#                                     SSI.input$min.Area,
#                                     SWSI.input$min.Area,
#                                     Ten.input$min.Area,
#                                     WR.input$min.Area,
#                                     Wool.input$min.Area), 
#                                 max(BR.input$max.Area,
#                                     CI.input$max.Area,
#                                     FR.input$max.Area,
#                                     FLI.input$max.Area,
#                                     IGS.input$max.Area,
#                                     JRFT.input$max.Area,
#                                     JRN.input$max.Area,
#                                     LEI.input$max.Area,
#                                     LMI.input$max.Area,
#                                     LL.input$max.Area,
#                                     Mud.input$max.Area,
#                                     NSI.input$max.Area,
#                                     NWSI.input$max.Area,
#                                     SSI.input$max.Area,
#                                     SWSI.input$max.Area,
#                                     Ten.input$max.Area,
#                                     WR.input$max.Area,
#                                     Wool.input$max.Area))) +
#   geom_line(aes(x = NWSIx.PLB, y = NWSIy.PLB), col = 'black', lwd = 1) +
#   annotate("text", x=10, y=10, label="North West Solitary Island") +
#   annotate("text", x=5, y=3, label = expression(paste(italic("b = "), -0.978))) +
#   theme_classic()
# 
# # Values of b to test to obtain confidence interval. For the real movement data
# # sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# # symmetric interval here.
# bvec = seq(PLB.bMLE.NWSI.b - 0.5, PLB.bMLE.NWSI.b + 0.5, 0.00001) 
# PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
# for(i in 1:length(bvec)){
#   PLB.LLvals[i] = negLL.PLB(bvec[i], x=NWSI.input$Area, n=length(NWSI.input$Area), xmin=NWSI.input$min.Area,
#                             xmax=NWSI.input$max.Area, sumlogx=NWSI.input$sum.log.Area)   
# }
# critVal = PLB.minNegLL.NWSI.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
# bIn95 = bvec[ PLB.LLvals < critVal ]
# NWSIbIn95 <- c(min(bIn95), max(bIn95))
# 
# # MLE South Solitary island (SSI) Area ===============================================================================
# # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# # as a starting point for nlm for MLE of b for PLB model.
# PLB.return.SSI <- mle_b(Site== "South Solitary Island", x=SSI.input$Area, log_x=SSI.input$log.Area, sum_log_x=SSI.input$sum.log.Area,
#                         x_min=SSI.input$min.Area, x_max=SSI.input$max.Area)
# PLB.bMLE.SSI.b <- PLB.return.SSI[[1]] 
# PLB.minLL.SSI.b <- PLB.return.SSI[[2]]
# 
# # plot and find 95% confidence intervals for MLE method.
# PLB.minNegLL.SSI.b <- PLB.minLL.SSI.b$minimum
# x <- SSI.input$Area
# SSIx.PLB = seq(min(SSI.input$Area), max(SSI.input$Area), length=1000) # x values to plot PLB. Note
# # that these encompass the data, and are not based
# # on the binning (in MEE Figure 6 the line starts as
# # min(x), not the first bin.
# SSIy.PLB = (1 - pPLB(x = SSIx.PLB, b = PLB.bMLE.SSI.b, xmin = min(SSIx.PLB),
#                      xmax = max(SSIx.PLB))) * length(SSI.input$Area)
# spectra.text <- as.character(round(PLB.bMLE.SSI.b, 2))
# SSIb_plot <- ggplot() +
#   geom_point(aes(x = (sort(SSI.input$Area, decreasing=TRUE)), y = (1:length(SSI.input$Area))),
#              color = "cadetblue", size = 2, alpha = 0.3) +
#   xlab(expression(paste("Colony area, ", italic("x"), " (cm^2)"))) +
#   ylab(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    "))) +
#   scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
#                      limits = c(0.25, max(length(BR.input$Area),
#                                           length(CI.input$Area),
#                                           length(FR.input$Area),
#                                           length(FLI.input$Area),
#                                           length(IGS.input$Area),
#                                           length(JRFT.input$Area),
#                                           length(JRN.input$Area),
#                                           length(LEI.input$Area),
#                                           length(LMI.input$Area),
#                                           length(LL.input$Area),
#                                           length(Mud.input$Area),
#                                           length(NSI.input$Area),
#                                           length(NWSI.input$Area),
#                                           length(SSI.input$Area),
#                                           length(SWSI.input$Area),
#                                           length(Ten.input$Area),
#                                           length(WR.input$Area),
#                                           length(Wool.input$Area)))) +
#   scale_x_continuous(trans = 'log10', breaks = c(0,1,5,10,100,1000,10000),
#                      limits = c(min(BR.input$min.Area,
#                                     CI.input$min.Area,
#                                     FR.input$min.Area,
#                                     FLI.input$min.Area,
#                                     IGS.input$min.Area,
#                                     JRFT.input$min.Area,
#                                     JRN.input$min.Area,
#                                     LEI.input$min.Area,
#                                     LMI.input$min.Area,
#                                     LL.input$min.Area,
#                                     Mud.input$min.Area,
#                                     NSI.input$min.Area,
#                                     NWSI.input$min.Area,
#                                     SSI.input$min.Area,
#                                     SWSI.input$min.Area,
#                                     Ten.input$min.Area,
#                                     WR.input$min.Area,
#                                     Wool.input$min.Area), 
#                                 max(BR.input$max.Area,
#                                     CI.input$max.Area,
#                                     FR.input$max.Area,
#                                     FLI.input$max.Area,
#                                     IGS.input$max.Area,
#                                     JRFT.input$max.Area,
#                                     JRN.input$max.Area,
#                                     LEI.input$max.Area,
#                                     LMI.input$max.Area,
#                                     LL.input$max.Area,
#                                     Mud.input$max.Area,
#                                     NSI.input$max.Area,
#                                     NWSI.input$max.Area,
#                                     SSI.input$max.Area,
#                                     SWSI.input$max.Area,
#                                     Ten.input$max.Area,
#                                     WR.input$max.Area,
#                                     Wool.input$max.Area))) +
#   geom_line(aes(x = SSIx.PLB, y = SSIy.PLB), col = 'black', lwd = 1) +
#   annotate("text", x=5, y=10, label="South Solitary Island") +
#   annotate("text", x=5, y=3, label = expression(paste(italic("b = "), -1.06))) +
#   theme_classic()
# 
# # Values of b to test to obtain confidence interval. For the real movement data
# # sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# # symmetric interval here.
# bvec = seq(PLB.bMLE.SSI.b - 0.5, PLB.bMLE.SSI.b + 0.5, 0.00001) 
# PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
# for(i in 1:length(bvec)){
#   PLB.LLvals[i] = negLL.PLB(bvec[i], x=SSI.input$Area, n=length(SSI.input$Area), xmin=SSI.input$min.Area,
#                             xmax=SSI.input$max.Area, sumlogx=SSI.input$sum.log.Area)   
# }
# critVal = PLB.minNegLL.SSI.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
# bIn95 = bvec[ PLB.LLvals < critVal ]
# SSIbIn95 <- c(min(bIn95), max(bIn95))
# 
# 
# # MLE South West Solitary Island (SWSI) Area ===============================================================================
# # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# # as a starting point for nlm for MLE of b for PLB model.
# PLB.return.SWSI <- mle_b(Site== "South West Solitary Island", x=SWSI.input$Area, log_x=SWSI.input$log.Area, sum_log_x=SWSI.input$sum.log.Area,
#                          x_min=SWSI.input$min.Area, x_max=SWSI.input$max.Area)
# PLB.bMLE.SWSI.b <- PLB.return.SWSI[[1]] 
# PLB.minLL.SWSI.b <- PLB.return.SWSI[[2]]
# 
# # plot and find 95% confidence intervals for MLE method.
# PLB.minNegLL.SWSI.b <- PLB.minLL.SWSI.b$minimum
# x <- SWSI.input$Area
# SWSIx.PLB = seq(min(SWSI.input$Area), max(SWSI.input$Area), length=1000) # x values to plot PLB. Note
# # that these encompass the data, and are not based
# # on the binning (in MEE Figure 6 the line starts as
# # min(x), not the first bin.
# SWSIy.PLB = (1 - pPLB(x = SWSIx.PLB, b = PLB.bMLE.SWSI.b, xmin = min(SWSIx.PLB),
#                       xmax = max(SWSIx.PLB))) * length(SWSI.input$Area)
# spectra.text <- as.character(round(PLB.bMLE.SWSI.b, 2))
# SWSIb_plot <- ggplot() +
#   geom_point(aes(x = (sort(SWSI.input$Area, decreasing=TRUE)), y = (1:length(SWSI.input$Area))),
#              color = "cadetblue", size = 2, alpha = 0.3) +
#   xlab(expression(paste("Colony area, ", italic("x"), " (cm^2)"))) +
#   ylab(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    "))) +
#   scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
#                      limits = c(0.25, max(length(BR.input$Area),
#                                           length(CI.input$Area),
#                                           length(FR.input$Area),
#                                           length(FLI.input$Area),
#                                           length(IGS.input$Area),
#                                           length(JRFT.input$Area),
#                                           length(JRN.input$Area),
#                                           length(LEI.input$Area),
#                                           length(LMI.input$Area),
#                                           length(LL.input$Area),
#                                           length(Mud.input$Area),
#                                           length(NSI.input$Area),
#                                           length(NWSI.input$Area),
#                                           length(SSI.input$Area),
#                                           length(SWSI.input$Area),
#                                           length(Ten.input$Area),
#                                           length(WR.input$Area),
#                                           length(Wool.input$Area)))) +
#   scale_x_continuous(trans = 'log10', breaks = c(0,1,5,10,100,1000,10000),
#                      limits = c(min(BR.input$min.Area,
#                                     CI.input$min.Area,
#                                     FR.input$min.Area,
#                                     FLI.input$min.Area,
#                                     IGS.input$min.Area,
#                                     JRFT.input$min.Area,
#                                     JRN.input$min.Area,
#                                     LEI.input$min.Area,
#                                     LMI.input$min.Area,
#                                     LL.input$min.Area,
#                                     Mud.input$min.Area,
#                                     NSI.input$min.Area,
#                                     NWSI.input$min.Area,
#                                     SSI.input$min.Area,
#                                     SWSI.input$min.Area,
#                                     Ten.input$min.Area,
#                                     WR.input$min.Area,
#                                     Wool.input$min.Area), 
#                                 max(BR.input$max.Area,
#                                     CI.input$max.Area,
#                                     FR.input$max.Area,
#                                     FLI.input$max.Area,
#                                     IGS.input$max.Area,
#                                     JRFT.input$max.Area,
#                                     JRN.input$max.Area,
#                                     LEI.input$max.Area,
#                                     LMI.input$max.Area,
#                                     LL.input$max.Area,
#                                     Mud.input$max.Area,
#                                     NSI.input$max.Area,
#                                     NWSI.input$max.Area,
#                                     SSI.input$max.Area,
#                                     SWSI.input$max.Area,
#                                     Ten.input$max.Area,
#                                     WR.input$max.Area,
#                                     Wool.input$max.Area))) +
#   geom_line(aes(x = SWSIx.PLB, y = SWSIy.PLB), col = 'black', lwd = 1) +
#   annotate("text", x=10, y=10, label="South West Solitary Island") +
#   annotate("text", x=5, y=3, label = expression(paste(italic("b = "), -0.967))) +
#   theme_classic()
# 
# # Values of b to test to obtain confidence interval. For the real movement data
# # sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# # symmetric interval here.
# bvec = seq(PLB.bMLE.SWSI.b - 0.5, PLB.bMLE.SWSI.b + 0.5, 0.00001) 
# PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
# for(i in 1:length(bvec)){
#   PLB.LLvals[i] = negLL.PLB(bvec[i], x=SWSI.input$Area, n=length(SWSI.input$Area), xmin=SWSI.input$min.Area,
#                             xmax=SWSI.input$max.Area, sumlogx=SWSI.input$sum.log.Area)   
# }
# critVal = PLB.minNegLL.SWSI.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
# bIn95 = bvec[ PLB.LLvals < critVal ]
# SWSIbIn95 <- c(min(bIn95), max(bIn95))
# 
# # MLE Tenenments (Ten) Area ===============================================================================
# # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# # as a starting point for nlm for MLE of b for PLB model.
# PLB.return.Ten <- mle_b(Site== "Tenements 1", x=Ten.input$Area, log_x=Ten.input$log.Area, sum_log_x=Ten.input$sum.log.Area,
#                         x_min=Ten.input$min.Area, x_max=Ten.input$max.Area)
# PLB.bMLE.Ten.b <- PLB.return.Ten[[1]] 
# PLB.minLL.Ten.b <- PLB.return.Ten[[2]]
# 
# # plot and find 95% confidence intervals for MLE method.
# PLB.minNegLL.Ten.b <- PLB.minLL.Ten.b$minimum
# x <- Ten.input$Area
# Tenx.PLB = seq(min(Ten.input$Area), max(Ten.input$Area), length=1000) # x values to plot PLB. Note
# # that these encompass the data, and are not based
# # on the binning (in MEE Figure 6 the line starts as
# # min(x), not the first bin.
# Teny.PLB = (1 - pPLB(x = Tenx.PLB, b = PLB.bMLE.Ten.b, xmin = min(Tenx.PLB),
#                      xmax = max(Tenx.PLB))) * length(Ten.input$Area)
# spectra.text <- as.character(round(PLB.bMLE.Ten.b, 2))
# Tenb_plot <- ggplot() +
#   geom_point(aes(x = (sort(Ten.input$Area, decreasing=TRUE)), y = (1:length(Ten.input$Area))),
#              color = "deeppink3", size = 2, alpha = 0.3) +
#   xlab(expression(paste("Colony area, ", italic("x"), " (cm^2)"))) +
#   ylab(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    "))) +
#   scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
#                      limits = c(0.25, max(length(BR.input$Area),
#                                           length(CI.input$Area),
#                                           length(FR.input$Area),
#                                           length(FLI.input$Area),
#                                           length(IGS.input$Area),
#                                           length(JRFT.input$Area),
#                                           length(JRN.input$Area),
#                                           length(LEI.input$Area),
#                                           length(LMI.input$Area),
#                                           length(LL.input$Area),
#                                           length(Mud.input$Area),
#                                           length(NSI.input$Area),
#                                           length(NWSI.input$Area),
#                                           length(SSI.input$Area),
#                                           length(SWSI.input$Area),
#                                           length(Ten.input$Area),
#                                           length(WR.input$Area),
#                                           length(Wool.input$Area)))) +
#   scale_x_continuous(trans = 'log10', breaks = c(0,1,5,10,100,1000,10000),
#                      limits = c(min(BR.input$min.Area,
#                                     CI.input$min.Area,
#                                     FR.input$min.Area,
#                                     FLI.input$min.Area,
#                                     IGS.input$min.Area,
#                                     JRFT.input$min.Area,
#                                     JRN.input$min.Area,
#                                     LEI.input$min.Area,
#                                     LMI.input$min.Area,
#                                     LL.input$min.Area,
#                                     Mud.input$min.Area,
#                                     NSI.input$min.Area,
#                                     NWSI.input$min.Area,
#                                     SSI.input$min.Area,
#                                     SWSI.input$min.Area,
#                                     Ten.input$min.Area,
#                                     WR.input$min.Area,
#                                     Wool.input$min.Area), 
#                                 max(BR.input$max.Area,
#                                     CI.input$max.Area,
#                                     FR.input$max.Area,
#                                     FLI.input$max.Area,
#                                     IGS.input$max.Area,
#                                     JRFT.input$max.Area,
#                                     JRN.input$max.Area,
#                                     LEI.input$max.Area,
#                                     LMI.input$max.Area,
#                                     LL.input$max.Area,
#                                     Mud.input$max.Area,
#                                     NSI.input$max.Area,
#                                     NWSI.input$max.Area,
#                                     SSI.input$max.Area,
#                                     SWSI.input$max.Area,
#                                     Ten.input$max.Area,
#                                     WR.input$max.Area,
#                                     Wool.input$max.Area))) +
#   geom_line(aes(x = Tenx.PLB, y = Teny.PLB), col = 'black', lwd = 1) +
#   annotate("text", x=5, y=10, label="Tenements 1") +
#   annotate("text", x=5, y=3, label = expression(paste(italic("b = "), -1.14))) +
#   theme_classic()
# 
# # Values of b to test to obtain confidence interval. For the real movement data
# # sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# # symmetric interval here.
# bvec = seq(PLB.bMLE.Ten.b - 0.5, PLB.bMLE.Ten.b + 0.5, 0.00001) 
# PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
# for(i in 1:length(bvec)){
#   PLB.LLvals[i] = negLL.PLB(bvec[i], x=Ten.input$Area, n=length(Ten.input$Area), xmin=Ten.input$min.Area,
#                             xmax=Ten.input$max.Area, sumlogx=Ten.input$sum.log.Area)   
# }
# critVal = PLB.minNegLL.Ten.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
# bIn95 = bvec[ PLB.LLvals < critVal ]
# TenbIn95 <- c(min(bIn95), max(bIn95))
# 
# 
# # MLE Wolf Rock (WR) Area ===============================================================================
# # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# # as a starting point for nlm for MLE of b for PLB model.
# PLB.return.WR <- mle_b(Site== "Wolf Rock", x=WR.input$Area, log_x=WR.input$log.Area, sum_log_x=WR.input$sum.log.Area,
#                        x_min=WR.input$min.Area, x_max=WR.input$max.Area)
# PLB.bMLE.WR.b <- PLB.return.WR[[1]] 
# PLB.minLL.WR.b <- PLB.return.WR[[2]]
# 
# # plot and find 95% confidence intervals for MLE method.
# PLB.minNegLL.WR.b <- PLB.minLL.WR.b$minimum
# x <- WR.input$Area
# WRx.PLB = seq(min(WR.input$Area), max(WR.input$Area), length=1000) # x values to plot PLB. Note
# # that these encompass the data, and are not based
# # on the binning (in MEE Figure 6 the line starts as
# # min(x), not the first bin.
# WRy.PLB = (1 - pPLB(x = WRx.PLB, b = PLB.bMLE.WR.b, xmin = min(WRx.PLB),
#                     xmax = max(WRx.PLB))) * length(WR.input$Area)
# spectra.text <- as.character(round(PLB.bMLE.WR.b, 2))
# WRb_plot <- ggplot() +
#   geom_point(aes(x = (sort(WR.input$Area, decreasing=TRUE)), y = (1:length(WR.input$Area))),
#              color = "cadetblue", size = 2, alpha = 0.3) +
#   xlab(expression(paste("Colony area, ", italic("x"), " (cm^2)"))) +
#   ylab(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    "))) +
#   scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
#                      limits = c(0.25, max(length(BR.input$Area),
#                                           length(CI.input$Area),
#                                           length(FR.input$Area),
#                                           length(FLI.input$Area),
#                                           length(IGS.input$Area),
#                                           length(JRFT.input$Area),
#                                           length(JRN.input$Area),
#                                           length(LEI.input$Area),
#                                           length(LMI.input$Area),
#                                           length(LL.input$Area),
#                                           length(Mud.input$Area),
#                                           length(NSI.input$Area),
#                                           length(NWSI.input$Area),
#                                           length(SSI.input$Area),
#                                           length(SWSI.input$Area),
#                                           length(Ten.input$Area),
#                                           length(WR.input$Area),
#                                           length(Wool.input$Area)))) +
#   scale_x_continuous(trans = 'log10', breaks = c(0,1,5,10,100,1000,10000),
#                      limits = c(min(BR.input$min.Area,
#                                     CI.input$min.Area,
#                                     FR.input$min.Area,
#                                     FLI.input$min.Area,
#                                     IGS.input$min.Area,
#                                     JRFT.input$min.Area,
#                                     JRN.input$min.Area,
#                                     LEI.input$min.Area,
#                                     LMI.input$min.Area,
#                                     LL.input$min.Area,
#                                     Mud.input$min.Area,
#                                     NSI.input$min.Area,
#                                     NWSI.input$min.Area,
#                                     SSI.input$min.Area,
#                                     SWSI.input$min.Area,
#                                     Ten.input$min.Area,
#                                     WR.input$min.Area,
#                                     Wool.input$min.Area), 
#                                 max(BR.input$max.Area,
#                                     CI.input$max.Area,
#                                     FR.input$max.Area,
#                                     FLI.input$max.Area,
#                                     IGS.input$max.Area,
#                                     JRFT.input$max.Area,
#                                     JRN.input$max.Area,
#                                     LEI.input$max.Area,
#                                     LMI.input$max.Area,
#                                     LL.input$max.Area,
#                                     Mud.input$max.Area,
#                                     NSI.input$max.Area,
#                                     NWSI.input$max.Area,
#                                     SSI.input$max.Area,
#                                     SWSI.input$max.Area,
#                                     Ten.input$max.Area,
#                                     WR.input$max.Area,
#                                     Wool.input$max.Area))) +
#   geom_line(aes(x = WRx.PLB, y = WRy.PLB), col = 'black', lwd = 1) +
#   annotate("text", x=5, y=10, label="Wolf Rock") +
#   annotate("text", x=5, y=3, label = expression(paste(italic("b = "), -0.974))) +
#   theme_classic()
# 
# # Values of b to test to obtain confidence interval. For the real movement data
# # sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# # symmetric interval here.
# bvec = seq(PLB.bMLE.WR.b - 0.5, PLB.bMLE.WR.b + 0.5, 0.00001) 
# PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
# for(i in 1:length(bvec)){
#   PLB.LLvals[i] = negLL.PLB(bvec[i], x=WR.input$Area, n=length(WR.input$Area), xmin=WR.input$min.Area,
#                             xmax=WR.input$max.Area, sumlogx=WR.input$sum.log.Area)   
# }
# critVal = PLB.minNegLL.WR.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
# bIn95 = bvec[ PLB.LLvals < critVal ]
# WRbIn95 <- c(min(bIn95), max(bIn95))
# 
# 
# # MLE Woolgoolga Reef (Wool) Area ===============================================================================
# # Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2007)
# # as a starting point for nlm for MLE of b for PLB model.
# PLB.return.Wool <- mle_b(Site== "Woolgoolga Reef", x=Wool.input$Area, log_x=Wool.input$log.Area, sum_log_x=Wool.input$sum.log.Area,
#                          x_min=Wool.input$min.Area, x_max=Wool.input$max.Area)
# PLB.bMLE.Wool.b <- PLB.return.Wool[[1]] 
# PLB.minLL.Wool.b <- PLB.return.Wool[[2]]
# 
# # plot and find 95% confidence intervals for MLE method.
# PLB.minNegLL.Wool.b <- PLB.minLL.Wool.b$minimum
# x <- Wool.input$Area
# Woolx.PLB = seq(min(Wool.input$Area), max(Wool.input$Area), length=1000) # x values to plot PLB. Note
# # that these encompass the data, and are not based
# # on the binning (in MEE Figure 6 the line starts as
# # min(x), not the first bin.
# Wooly.PLB = (1 - pPLB(x = Woolx.PLB, b = PLB.bMLE.Wool.b, xmin = min(Woolx.PLB),
#                       xmax = max(Woolx.PLB))) * length(Wool.input$Area)
# spectra.text <- as.character(round(PLB.bMLE.Wool.b, 2))
# Woolb_plot <- ggplot() +
#   geom_point(aes(x = (sort(Wool.input$Area, decreasing=TRUE)), y = (1:length(Wool.input$Area))),
#              color = "darkblue", size = 2, alpha = 0.3) +
#   xlab(expression(paste("Colony area, ", italic("x"), " (cm^2)"))) +
#   ylab(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    "))) +
#   scale_y_continuous(trans = 'log10', breaks = c(1,10,100,500,3000),
#                      limits = c(0.25, max(length(BR.input$Area),
#                                           length(CI.input$Area),
#                                           length(FR.input$Area),
#                                           length(FLI.input$Area),
#                                           length(IGS.input$Area),
#                                           length(JRFT.input$Area),
#                                           length(JRN.input$Area),
#                                           length(LEI.input$Area),
#                                           length(LMI.input$Area),
#                                           length(LL.input$Area),
#                                           length(Mud.input$Area),
#                                           length(NSI.input$Area),
#                                           length(NWSI.input$Area),
#                                           length(SSI.input$Area),
#                                           length(SWSI.input$Area),
#                                           length(Ten.input$Area),
#                                           length(WR.input$Area),
#                                           length(Wool.input$Area)))) +
#   scale_x_continuous(trans = 'log10', breaks = c(0,1,5,10,100,1000,10000),
#                      limits = c(min(BR.input$min.Area,
#                                     CI.input$min.Area,
#                                     FR.input$min.Area,
#                                     FLI.input$min.Area,
#                                     IGS.input$min.Area,
#                                     JRFT.input$min.Area,
#                                     JRN.input$min.Area,
#                                     LEI.input$min.Area,
#                                     LMI.input$min.Area,
#                                     LL.input$min.Area,
#                                     Mud.input$min.Area,
#                                     NSI.input$min.Area,
#                                     NWSI.input$min.Area,
#                                     SSI.input$min.Area,
#                                     SWSI.input$min.Area,
#                                     Ten.input$min.Area,
#                                     WR.input$min.Area,
#                                     Wool.input$min.Area), 
#                                 max(BR.input$max.Area,
#                                     CI.input$max.Area,
#                                     FR.input$max.Area,
#                                     FLI.input$max.Area,
#                                     IGS.input$max.Area,
#                                     JRFT.input$max.Area,
#                                     JRN.input$max.Area,
#                                     LEI.input$max.Area,
#                                     LMI.input$max.Area,
#                                     LL.input$max.Area,
#                                     Mud.input$max.Area,
#                                     NSI.input$max.Area,
#                                     NWSI.input$max.Area,
#                                     SSI.input$max.Area,
#                                     SWSI.input$max.Area,
#                                     Ten.input$max.Area,
#                                     WR.input$max.Area,
#                                     Wool.input$max.Area))) +
#   geom_line(aes(x = Woolx.PLB, y = Wooly.PLB), col = 'black', lwd = 1) +
#   annotate("text", x=5, y=10, label="Woolgoolga Reef") +
#   annotate("text", x=5, y=3, label = expression(paste(italic("b = "), -0.806))) +
#   theme_classic()
# 
# # Values of b to test to obtain confidence interval. For the real movement data
# # sets in Table 2 of Edwards (2011) the intervals were symmetric, so make a
# # symmetric interval here.
# bvec = seq(PLB.bMLE.Wool.b - 0.5, PLB.bMLE.Wool.b + 0.5, 0.00001) 
# PLB.LLvals = vector(length=length(bvec))  # negative log-likelihood for bvec
# for(i in 1:length(bvec)){
#   PLB.LLvals[i] = negLL.PLB(bvec[i], x=Wool.input$Area, n=length(Wool.input$Area), xmin=Wool.input$min.Area,
#                             xmax=Wool.input$max.Area, sumlogx=Wool.input$sum.log.Area)   
# }
# critVal = PLB.minNegLL.Wool.b  + qchisq(0.95,1)/2 # 1 degree of freedom, Hilborn and Mangel (1997) p162.
# bIn95 = bvec[ PLB.LLvals < critVal ]
# WoolbIn95 <- c(min(bIn95), max(bIn95))
# 
# 
# # Plot all the size spectra slopes ==============================================
# all_SS <- cowplot::plot_grid(BRb_plot+theme(axis.title = element_blank()), SSIb_plot+theme(axis.title = element_blank()), 
#                              SWSIb_plot+theme(axis.title = element_blank()), Woolb_plot+theme(axis.title = element_blank()),
#                              NWSIb_plot+theme(axis.title = element_blank()), NRb_plot+theme(axis.title = element_blank()), NSIb_plot+theme(axis.title = element_blank()),
#                              JRFTb_plot+theme(axis.title = element_blank()), JRNb_plot+theme(axis.title = element_blank()), 
#                              CIb_plot+theme(axis.title = element_blank()), FRb_plot+theme(axis.title = element_blank()), Henb_plot+theme(axis.title = element_blank()), 
#                              FLIb_plot+theme(axis.title = element_blank()), IGSb_plot+theme(axis.title = element_blank()),
#                              Mudb_plot+theme(axis.title = element_blank()), WRb_plot+theme(axis.title = element_blank()), 
#                              LEIb_plot+theme(axis.title = element_blank()), LMIb_plot+theme(axis.title = element_blank()), 
#                              LLb_plot+theme(axis.title = element_blank()), Tenb_plot+theme(axis.title = element_blank()), 
#                              ncol = 4, nrow = 5, labels="AUTO",label_size = 12, hjust = c(0.5,0.5,0.5,
#                                                                                           -0.5,-0.5,-0.5,
#                                                                                           -0.5,-0.5,-0.5,
#                                                                                           -0.5,-0.5,-0.5,
#                                                                                           -0.5,-0.5,-0.5,
#                                                                                           -0.5,-0.5,-0.5,-0.5,-0.5), 
#                              vjust=c(1.5,1.5,1.5,1.5,
#                                      0.2,0.2,0.2,0.2,
#                                      0.2,0.2,0.2,0.2,
#                                      0.2,0.2,0.2,0.2,
#                                      0.2,0.2,0.2,0.2), align="hv")
# 
# #create common x and y labels https://stackoverflow.com/questions/33114380/centered-x-axis-label-for-muliplot-using-cowplot-package
# y.grob <- grid::textGrob(expression(paste("Number of colonies with sizes", " ">=" ", italic("x"), "    ")), rot=90)
# x.grob <- grid::textGrob(expression(paste("Colony area, ", italic("x"), (cm^"2"))))
# 
# #add to plot
# spectra <- gridExtra::grid.arrange(all_SS, left = y.grob, bottom = x.grob)
# 
# ggsave("spectra.svg", bg = "transparent", width = 350, height = 210, units = "mm")
