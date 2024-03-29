# set.params - create list with parameters for running MLE
#Code adopted from Edwards et al. (2017).
set.params = function(Area){
  # Creates a list of parameters for input into the negll.PLB and pPLB functions
  Area <- Area
  log.Area <- log(Area)
  sum.log.Area <- sum(log.Area)
  min.Area <- min(Area)
  max.Area <- max(Area)
  n <- length(Area)
  out.list <- (list(Area, log.Area, sum.log.Area, min.Area, max.Area, n))
  names(out.list) <- c("Area", "log.Area", "sum.log.Area", "min.Area", "max.Area", "n")
  return(out.list)
}

# Use analytical value of MLE b for PL model (Box 1, Edwards et al. 2017)
# as a starting point for nlm for MLE of b for PLB model. Code adopted from
# Edwards et al. (2017).
mle_b = function(region, x, log_x, sum_log_x, x_min, x_max){
  
  PL.bMLE = 1/( log(min(x)) - sum_log_x/length(x)) - 1
  
  PLB.minLL =  nlm(negLL.PLB, p=PL.bMLE, x=x, n=length(x),
                   xmin=x_min, xmax=x_max, sumlogx=sum_log_x, hessian = TRUE) #, print.level=2 )
  
  PLB.bMLE = PLB.minLL$estimate
  
  PLB.return = list(PLB.bMLE, PLB.minLL)
  
  return(PLB.return)
}

# negLL.PLB - negative log-likelihood function (function by Edwards et al. 2017)
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

# pPLB - bounded power-law distribution function (function by Edwards et al. 2017)
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

#calculate g(x), the probability we minus-sample a circle of area x
#Arguments:
#w, v: rectangular window dimensions
#x: circle area
#Value: probability a circle of area x is minus-sampled
#NOTE: IGNORES THE DENOMINATOR w * v, which cancels out in the formulas for minus-sampled distributions
getgx <- function(v, w, x){
  return(v * w - 2 * (v + w) * sqrt(x / pi) + 4 * x / pi)
}

#density for minus-sampled bounded power law (in form we can supply to optimizer)
#Arguments:
#b: power law coefficient
#x: value at which to get density
#xmin: min for bounded power law
#w, v: width and height of sampling window (ASSUMES v <= w and xmax greater than largest observable object)
#Value: density f(x) at x
dMSBPLanalytic <- function(b, x, xmin, w, v){
  numerator <- 4 / pi * x ^ (b + 1) - 2 * (w + v) / sqrt(pi) * x ^ (b + 0.5) + w * v * x ^ b
  xmw <- pi / 4 * v ^ 2 #largest circle we can fit in window
  denominator <- integrateMSBPL(x = xmw, b = b, xmin = xmin, w = w, v = v)
  fx <- numerator / denominator
  fx[x < xmin] <- 0
  fx[x > xmw] <- 0
  return(fx)
}

#CDF for minus sampled bounded power law
#Arguments:
#b: power law coefficient
#x: value at which to get density
#xmin: min for bounded power law
#w, v: width and height of sampling window (ASSUMES v <= w and xmax greater than largest observable object)
#Value: cdf for minus sampled bounded power law up to area x
FMSBPL <- function(b, x, xmin, w, v){
  xmw <- pi / 4 * v ^ 2 #largest circle we can fit in window
  denominator <- integrateMSBPL(x = xmw, b = b, xmin = xmin, w = w, v = v)
  numerator <- integrateMSBPL(x = x, b = b, xmin = xmin, w = w, v = v)
  return(numerator/denominator)
}

#inverse CDF for minus sampled bounded power law
# vector of u in [0,1]
#b: power law coefficient
#xmin: min for bounded power law
#w, v: width and height of sampling window (ASSUMES v <= w and xmax greater than largest observable object)
#Value: inverse CDF evaluated at u
FMSBPLinv <- function(u, b, xmin, w, v){
  nu <- length(u)
  xout <- numeric(nu)
  xmw <- pi / 4 * v ^ 2 
  for(i in 1:nu){
    xout[i] <- uniroot(f = function(x){
      FMSBPL(b = b, x = x, xmin = xmin, w = w, v = v) - u[i]
     }, lower = xmin, upper = xmw)$root
  }
  return(xout) #evaluated for each point in u
}

#integrate unstandardised minus sampled bounded power law from xmin to x
#Arguments:
#x: upper bound for integral
#b: power law coefficient
#xmin: min for bounded power law
#w, v: width and height of sampling window (ASSUMES v <= w and xmax greater than largest observable object)
#Value: the integral
integrateMSBPL <- function(x, b, xmin, w, v){
  if(b == -2){#special case for denominator
    dplus <- w * v / (b + 1) * x ^ (b + 1) - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * x ^ (b + 1.5) + 4 / pi * log(x)
    dminus <- w * v / (b + 1) * xmin ^ (b + 1) - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * xmin ^ (b + 1.5) + 4 / pi * log(xmin)
  } else if (b == - 1.5){#special case for denominator
    dplus <- w * v / (b + 1) * x ^ (b + 1) - 2 * (w + v) / sqrt(pi) * log(x) + 4 / (pi * (b + 2)) * x ^ (b + 2)
    dminus <- w * v / (b + 1) * xmin ^ (b + 1) - 2 * (w + v) / sqrt(pi) * log(xmin) + 4 / (pi * (b + 2)) * xmin ^ (b + 2)
  } else if (b == -1){#special case for denominator
    dplus <- w * v * log(x) - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * x ^ (b + 1.5) + 4 / (pi * (b + 2)) * x ^ (b + 2)
    dminus <- w * v * log(xmin) - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * xmin ^ (b + 1.5) + 4 / (pi * (b + 2)) * xmin ^ (b + 2)
  } else {
    dplus <- w * v / (b + 1) * x ^ (b + 1) - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * x ^ (b + 1.5) + 4 / (pi * (b + 2)) * x ^ (b + 2)
    dminus <- w * v / (b + 1) * xmin ^ (b + 1) - 2 * (w + v) / (sqrt(pi) * (b + 1.5)) * xmin ^ (b + 1.5) + 4 / (pi * (b + 2)) * xmin ^ (b + 2)
  }
  return(dplus - dminus)
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
dMSBPL <- function(x, b, C, xmin, xmax, w, v){
  numerator <- 4 / pi * x ^ (b + 1) - 2 * (w + v) / sqrt(pi) * x ^ (b + 0.5) + w * v * x ^ b
  n2 <- dboundedpowerlaw(x = x, b = b, C = C, xmin = xmin, xmax = xmax) * getpir(w = w, v = v, r = sqrt(x / pi))
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
#Argument: area x FC: should be in cm^2
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

#log likelihood and AIC for a lognormal distribution
#Arguments: x: vector of sizes (not logged)
#Value: list containing lognormal log likelihood and AIC
lnormAIC <- function(x){
  theta <- estimatelognormal(x = x)
  lllognorm <- sum(dlnorm(x = x, meanlog = theta$meanlog, sdlog = theta$sdlog, log = TRUE))
  AIC <- 2*2 - 2*lllognorm
  return(list(lllognorm = lllognorm , AIClognorm = AIC)) 
}

#maximum likelihood estimates of parameters of lognormal distribution
#Arguments: vector x of positive observations
#Value: list containing meanlog (ML estimate of mean of log x) and sdlog (ML estimate of sd of log x)
estimatelognormal <- function(x){
  n <- length(x)
  logx <- log(x)
  meanlog <- mean(logx)
  sdlog <- sd(logx) * sqrt((n - 1) / n) #ML estimate
  return(list(meanlog = meanlog, sdlog = sdlog))
}
# log likelihood and AIC of a bounded power law distribution ---- , AIC is probably wrong because log likelihood not locally quadratic around xmin

BPLAIC <- function(C, b, x){#making the argument be x instead of a may be easier to remember (same as for normal)
  n <- length(x)
  llBPL <- n * log(C) + b * sum(log(x)) 
  AIC <- 2*3 - 2*llBPL
  return(list(llBPL = llBPL, AICBPL = AIC))
}

#density for minus-sampled lognormal
#Arguments:
#x: value up to which we want to integrate
#mu, sigma: mean and sd of log area
#w, v: width and height of sampling window (ASSUMES v < w)
#log (default FALSE): return log density?
#Value: minus-sampled lognormal density
dMSlnorm <- function(x, mu, sigma, v, w, log = FALSE){
  xmw <- pi / 4 * v ^ 2 #largest circle we can fit in window
  if(log){
    lnum <- dlnorm(x = x, meanlog = mu, sdlog = sigma, log = TRUE) + log(getgx(v = v, w = w, x = x))
    ldenom <- log(dMSlnormintegral(x = xmw, mu = mu, sigma = sigma, w = w, v = v)$value)
    return(lnum - ldenom)
  } else {
    numerator <- dlnorm(x = x, meanlog = mu, sdlog = sigma) * getgx(v = v, w = w, x = x)
    denominator <- dMSlnormintegral(x = xmw, mu = mu, sigma = sigma, w = w, v = v)$value
    return(numerator / denominator)
  }
}

#Denominator for minus-sampled lognormal by numerical integration
#Arguments:
#x: value up to which we want to integrate
#mu, sigma: mean and sd of log area
#w, v: width and height of sampling window (ASSUMES v < w)
#Value: numerical integral of minus-sampled lognormal density from 0 to x
dMSlnormintegral <- function(x, mu, sigma, w, v){
  myfunction <- function(x){
    f <- dlnorm(x = x, meanlog = mu, sdlog = sigma) * getgx(w = w, v = v, x = x)
    return(f)
  }
  return(integrate(f = myfunction, lower = 0, upper = x))
}

#cdf for minus sampled lognormal
#Arguments:
#x: value at which to get cdf
#mu, sigma: mean and sd for log area
#w, v: width and height of sampling window (ASSUMES v <= w and xmax greater than largest observable object)
#Value: cdf for minus sampled lognormal up to area x
FMSlnorm <- function(x, mu, sigma, w, v){
  xmw <- pi / 4 * v ^ 2 #largest circle we can fit in window
  denominator <- dMSlnormintegral(x = xmw, mu = mu, sigma = sigma, w = w, v = v)$value
  numerator <- dMSlnormintegral(x = x, mu = mu, sigma = sigma, w = w, v = v)$value
  return(numerator / denominator)
}

#inverse CDF for minus sampled lognormal
#u in [0,1]
#mu, sigma: mean and sd for log area
#w, v: width and height of sampling window (ASSUMES v <= w and xmax greater than largest observable object)
#Value: inverse CDF evaluated at u
#Note: not vectorized. Use sapply() if we have a vector of values of u
FMSlnorminv <- function(u, mu, sigma, w, v){
  nu <- length(u)
  xout <- numeric(nu)
  xmw <- pi / 4 * v ^ 2 
  for(i in 1:nu){
    xout[i] <- uniroot(f = function(x){
      FMSlnorm(x = x, mu = mu, sigma = sigma, w = w, v = v) - u[i]
    }, lower = 0, upper = xmw)$root
  }
  return(xout)
}

#negative log likelihood for minus-sampled lognormal (in form we can supply to optimizer)
#Arguments:
#theta: parameter vector (mu, sigma): mean and sd of log area
#x: vector of sizes
#w, v: width and height of sampling window (ASSUMES v <= w and xmax greater than largest observable object)
#Value: negative log likelihood for observations x, with parameters mu, sigma
negloglikMSlnorm <- function(theta, x, w, v){
  logfx <- sum(dMSlnorm(x = x, mu = theta[1], sigma = theta[2], w = w, v = v, log = TRUE))
  return(-logfx)
}

#plot contours of negative log likelihood for minus-sampled lognormal
#Arguments:
#x: vector of sizes
#w, v: width and height of sampling window (ASSUMES v <= w and xmax greater than largest observable object)
#mu, sigma: true mean and sd of log size
#muhat, thetahat: estimated parameters
#Value: contours of negative log likelihood, with dot at true parameter vector
plotnegloglikMSlnorm <- function(x, w, v, mu, sigma, muhat, sigmahat){
  mycolors <- brewer.pal(3, "Dark2")
    muv <- seq(from = mu - 2, to = mu + 2, length.out = 40)
  sigmav <- seq(from = sigma / 2, to = 2 * sigma, length.out = 40)
  nll <- array(dim = c(length(muv), length(sigmav)))
  for(i in 1:length(muv)){
    for(j in 1:length(sigmav))
      nll[i, j] <- negloglikMSlnorm(theta = c(muv[i], sigmav[j]), x = x, w = w, v = v)
  }
   contour(muv, sigmav, nll, xlab = expression(mu), ylab = expression(sigma))
   points(mu, sigma, pch = 16, col = adjustcolor(mycolors[1], 0.4))
   points(muhat, sigmahat, pch = 1, col = adjustcolor(mycolors[2], 0.4))
   legend("topright", pch = c(16, 1), col = mycolors[1:2], legend = c("true", "estimated"), bty = "n")
}

#maximum likelihood estimate of parameters for minus-sampled lognormal
#x: vector of sizes
#w: width of window
#v: height of window
#Value: object returned by optim(). Contains par (parameter vector: mean and sd of log sizes), value (negative log likelihood), convergence (0 indicates success)
estimateMSlnorm <- function(x, w, v){
  testfun <- function(theta, x){#test: ordinary lognormal
    -sum(dlnorm(meanlog = theta[1], sdlog = theta[2], x = x, log = TRUE))
  }
  
  par <- c(mean(log(x)), sd(log(x))) #plausible initial guesses: sample mean and sd of log sizes
  #print("Ordinary lognormal:")
  #ordinary <- optim(f = testfun, par = par, method = "Nelder-Mead", x = x) #test: ordinary lognormal
  #print(c(ordinary$par, ordinary$value, ordinary$convergence))
  thetamslnorm <- optim(f = negloglikMSlnorm, par = par, method = "Nelder-Mead", x = x, w = w, v = v) #BFGS didn't behave nicely in this case
  return(thetamslnorm)
}

# AIC for minus sampled lognormal
# AIC = 2K - 2ln(L)
# K = number of model parameters - we have 2 for normal distribution
# ln(L) is the log likelihood from thetaML$value is the object returned by estimateMSlnorm

MSlnormAIC <- function(thetaML){
  AIC <- 2*2 - 2*(-thetaML$value)
  return(list(llmslognorm = -thetaML$value , AICmslognorm = AIC))
}

#AIC for minus sampled bounded power law, this is probably wrong because log likelihood not locally quadratic around xmin 
#arguments: msbplfit is an object returned by estimatebMSBPL()
MSBPLAIC <- function(msbplfit){
  AIC <- 2*2 - 2*(-msbplfit$objective)
  return(list(llMSBPL =-msbplfit$objective, AICMSBPL = AIC))
}
