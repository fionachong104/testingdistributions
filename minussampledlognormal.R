rm(list = ls())
library(RColorBrewer)

source("coralsizedistfuncs.R")

# images dimensions from CANON camera (w x v = 3648 x 2736)
# 10 cm = roughly 350 pixels; 1 m = 3500; 1 cm = 35; so dividing by 35 to see what that translates to in cm x cm
w <- 3648/35
v <- 2736/35

xmaxminus <- pi / 4 * v^2 #max area of circle that can fit in window
n <- 1e4 #number of circles
nplot <- 100 #number to plot
mu <- 5 #mean log area
sigma <- 1 #sd log area
x <- rlnorm(n = n, meanlog = mu, sdlog = sigma) #simulated areas from lognormal

r <- radiusfromarea(x) #colony radii (assume circles)
mycolors <- brewer.pal(3, "Dark2")

alpha <- runif(n = n, min = 0, max = w) #(alpha, beta) uniform random points in sampling window
beta <- runif(n = n, min = 0, max = v)

par(mfrow = c(3, 3))
plot(c(0, w), c(0, v), type = "n", asp = 1, xlab = "", ylab = "")
plotframe(w = w, v = v)
isinframe <- inframe(alpha = alpha, beta = beta, r = r, w = w, v = v)
for(i in 1:nplot){#drawing all the circles can be slow so make nplot not too large
  plotcircle(alpha = alpha[i], beta = beta[i], r = r[i], isinframe = isinframe[i])
}

hist(x, freq = FALSE, main = "lognormal", breaks = 100, xlab = expression(italic(x)), ylab = "density")
xseq <- seq(from = min(x), to = max(x), length.out = 1e4)
fx <- dlnorm(x = xseq, meanlog = mu, sdlog = sigma)
lines(xseq, fx, col = mycolors[1], lwd = 2, lty = "solid") #density curve for lognormal
legend("topright", bty = "n", lwd = 2, lty = c("solid"), col = c(mycolors[1]), legend = c("calculated"))

hist(x[isinframe], freq = FALSE, main = "minus-sampled lognormal", breaks = 100, xlab = expression(italic(x)), ylab = "density")
xseqMS <- seq(from = min(x), to = max(x), length.out = 1e4)
fxMS <- dMSlnorm(x = xseqMS, mu = mu, sigma = sigma, v = v, w = w)
lines(xseqMS, fxMS, col = mycolors[2], lwd = 2, lty = "dashed") #density curve for minus-sampled lognormal
legend("topright", bty = "n", lwd = 2, lty = c("dashed"), col = c(mycolors[2]), legend = c("calculated"))

plot(log(xseq), log(fx), type = "l", lwd = 2, lty = "solid", col = mycolors[1], xlab = expression(log(italic(x))), ylab = "log density")
lines(log(xseqMS), log(fxMS), col = mycolors[2], lwd = 2, lty = "dashed")
legend("bottomleft", bty = "n", lwd = 2, lty = c("solid", "dashed"), col = c(mycolors[1], mycolors[2]), legend = c("lognormal", "minus-sampled lognormal"))

#plotnegloglik(x = x[isinframe], w = w, v = v)
#abline(v = b, lty = "solid", col = mycolors[1])
#bML <- estimatebMSBPL(x = x[isinframe], w = w, v = v)
#abline(v = bML$minimum, lty = "dashed", col = mycolors[2])
#legend("topright", lty = c("solid", "dashed"), legend = c("true", "ML"), col = mycolors[1:2], bty = "n")

ninframe <- length(x[isinframe])
#sitex.MSBPL = seq(min(x[isinframe]), max(x[isinframe]), length = 1000)
#sitey.MSBPL = (1 - FMSBPL(x = sitex.MSBPL, b = b, xmin = xmin,
#                      w = w, v = v)) * ninframe
#plot(sort(x[isinframe], decreasing=TRUE),1:ninframe, xlab = expression(italic(x)), ylab = expression(italic(S(x)*n)[inframe]))
#lines(sitex.MSBPL,sitey.MSBPL)

#minus-sampled lognormal qqplot
qqplot(FMSlnorminv(u = ppoints(ninframe), mu = mu, sigma = sigma, w = w, v = v), x[isinframe], xlab = "theoretical quantiles", ylab = "sample quantiles", main = "minus-sampled lognormal Q-Q plot")

qqline(x[isinframe], distribution = function(p){
  FMSlnorminv(p, mu = mu, sigma = sigma, w = w, v = v)
})
