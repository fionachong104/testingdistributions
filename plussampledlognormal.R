rm(list = ls())
library(RColorBrewer)
#set.seed(123456789) #gives corner case check with rconst <- 10, n <- 10

source("coralsizedistfuncs.R")

# images dimensions from CANON camera (w x v = 3648 x 2736)
# 10 cm = roughly 350 pixels; 1 m = 3500; 1 cm = 35; so dividing by 35 to see what that translates to in cm x cm
w <- 3648/35
v <- 2736/35

xmaxminus <- pi / 4 * v^2 #max area of circle that can fit in window
maxdist <- 50 #max perpendicular distance from edge of image frame to simulate
n <- 1e3
rconst <- 10
r <- rep(rconst, n) #fixed r for now
nplot <- 1e2
mycolors <- brewer.pal(3, "Dark2")
alpha <- runif(n = n, min = 0 - maxdist, max = w + maxdist) #(alpha, beta) uniform random points in plus-sampling window
beta <- runif(n = n, min = 0 - maxdist, max = v + maxdist)
par(mfrow = c(1, 1))
plot(c(0 - maxdist, w + maxdist), c(0 - maxdist, v + maxdist), type = "n", asp = 1, xlab = "", ylab = "")
plotframe(w = w, v = v)
plotplusframe(w = w, v = v, r = rconst)
istrunc <- logical(n)
for(i in 1:n){
  istrunc[i] <- intrunc(alpha = alpha[i], beta = beta[i], r = r[i], w = w, v = v) #partly visible in frame NEEDS VECTORIZING
}
isinframe <- inframe(alpha = alpha, beta = beta, r = r, w = w, v = v) #entirely in frame
isvisible <- istrunc | isinframe #partly or entirely visible
for(i in 1:nplot){#drawing all the circles can be slow so make nplot not too large
  plotcircle(alpha = alpha[i], beta = beta[i], r = r[i], isinframe = istrunc[i])
}

#is probability of being truncated correct?
print(paste("simulated probability:", sum(istrunc) / sum(isvisible)))
if(v > 2 * rconst){
  P <- v * w + 2 * rconst * (v + w) + pi * rconst^2
  M <- (w - 2 * rconst) * (v - 2 * rconst)
  expectprob <- (P - M) / P
} else {
  expectprob <- 1
}
print(paste("expected probability:", expectprob))
      

#n <- 1e4 #number of circles
#nplot <- 100 #number to plot
#mu <- 5 #mean log area
#sigma <- 1 #sd log area
#x <- rlnorm(n = n, meanlog = mu, sdlog = sigma) #simulated areas from lognormal

#r <- radiusfromarea(x) #colony radii (assume circles)

#par(mfrow = c(3, 3))
#plot(c(0, w), c(0, v), type = "n", asp = 1, xlab = "", ylab = "")
#plotframe(w = w, v = v)
#isinframe <- inframe(alpha = alpha, beta = beta, r = r, w = w, v = v)
#for(i in 1:nplot){#drawing all the circles can be slow so make nplot not too large
#  plotcircle(alpha = alpha[i], beta = beta[i], r = r[i], isinframe = isinframe[i])
#}

#hist(x, freq = FALSE, main = "lognormal", breaks = 100, xlab = expression(italic(x)), ylab = "density")
#xseq <- seq(from = min(x), to = max(x), length.out = 1e4)
#fx <- dlnorm(x = xseq, meanlog = mu, sdlog = sigma)
#lines(xseq, fx, col = mycolors[1], lwd = 2, lty = "solid") #density curve for lognormal
#legend("topright", bty = "n", lwd = 2, lty = c("solid"), col = c(mycolors[1]), legend = c("calculated"))

#hist(x[isinframe], freq = FALSE, main = "minus-sampled lognormal", breaks = 100, xlab = expression(italic(x)), ylab = "density")
#xseqMS <- seq(from = min(x), to = max(x), length.out = 1e4)
#fxMS <- dMSlnorm(x = xseqMS, mu = mu, sigma = sigma, v = v, w = w)
#lines(xseqMS, fxMS, col = mycolors[2], lwd = 2, lty = "dashed") #density curve for minus-sampled lognormal
#legend("topright", bty = "n", lwd = 2, lty = c("dashed"), col = c(mycolors[2]), legend = c("calculated"))

#plot(log(xseq), log(fx), type = "l", lwd = 2, lty = "solid", col = mycolors[1], xlab = expression(log(italic(x))), ylab = "log density")
#lines(log(xseqMS), log(fxMS), col = mycolors[2], lwd = 2, lty = "dashed")
#legend("bottomleft", bty = "n", lwd = 2, lty = c("solid", "dashed"), col = c(mycolors[1], mycolors[2]), legend = c("lognormal", "minus-sampled lognormal"))

#thetaML <- estimateMSlnorm(x = x[isinframe], w = w, v = v)
#plotnegloglikMSlnorm(x = x[isinframe], w = w, v = v, mu = mu, sigma = sigma, muhat = thetaML$par[1], sigmahat = thetaML$par[2]) #contours of negative log likelihood
#print("Parameter estimates:")
#print(thetaML$par)

#ninframe <- length(x[isinframe])
#sitex.MSlnorm <- seq(min(x[isinframe]), max(x[isinframe]), length = 1000)
#sitey.MSlnorm <- (1 - sapply(sitex.MSlnorm, FUN = FMSlnorm, mu = mu, sigma = sigma, w = w, v = v)) * ninframe
#plot(sort(x[isinframe], decreasing=TRUE), 1:ninframe, xlab = expression(italic(x)), ylab = expression(italic(S(x)*n)[inframe]))
#lines(sitex.MSlnorm, sitey.MSlnorm)

#minus-sampled lognormal qqplot
#qqplot(FMSlnorminv(u = ppoints(ninframe), mu = mu, sigma = sigma, w = w, v = v), x[isinframe], xlab = "theoretical quantiles", ylab = "sample quantiles", main = "minus-sampled lognormal Q-Q plot")

#qqline(x[isinframe], distribution = function(p){
#  FMSlnorminv(p, mu = mu, sigma = sigma, w = w, v = v)
#})
