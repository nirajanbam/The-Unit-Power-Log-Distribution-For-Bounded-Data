
library(goftest)
library(maxLik)
y <- c(
  0.5935, 0.6623, 0.2917, 0.7705, 0.5000, 0.4839, 0.7536, 0.9048, 0.8667, 0.9189,
  0.9118, 0.7660, 0.5111, 0.2162, 0.8485, 0.8471, 0.9091, 0.2308, 0.7209, 0.9000,
  0.8000, 0.4828, 0.3333, 0.8605, 0.1250, 1.0000, 0.0000, 0.0000, 0.9870, 1.0000,
  1.0000, 0.5000, 0.9167, 1.0000, 1.0000, 1.0000, 1.0000, 0.5714, 0.9826, 0.8571,
  1.0000, 1.0000, 0.9714, 0.8889, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.3000,
  0.9344, 1.0000, 0.9831, 1.0000, 0.9565, 0.9833, 0.9500, 0.0000, 1.0000, 1.0000,
  1.0000, 1.0000, 1.0000, 0.9630, 1.0000, 1.0000, 1.0000, 1.0000, 0.9231, 0.6061,
  0.9048, 1.0000, 0.9667, 0.9583, 0.0000, 1.0000, 1.0000, 0.5636, 0.9878, 0.2000,
  0.9600, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.7419, 0.9524, 1.0000, 1.0000,
  0.8182, 0.5000, 0.9545, 0.9375, 0.8750, 0.8929, 1.0000, 1.0000, 0.9714, 1.0000,
  1.0000, 1.0000, 1.0000, 0.8000, 1.0000, 0.6500, 0.6774, 0.9412, 0.0556, 0.5106,
  1.0000, 0.1948, 1.0000, 1.0000, 0.8824, 1.0000, 1.0000, 1.0588, 0.8571, 1.0435,
  0.5833, 1.0000, 0.6875, 0.8182, 1.0000, 1.0000, 0.8667, 1.0000, 0.9167, 0.9000,
  1.0000, 1.0000, 1.0000, 0.7778, 1.0000, 0.9000, 0.9231, 1.0000, 0.8750, 0.9286,
  1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 0.5000, 0.8750, 1.0000, 0.8929, 0.5000,
  1.0000, 1.0000, 0.9545
)

x<-  y[1:29]        # y[1:29] # p-value 0.87
n<- length(x)
hist(x)
range(x)


################# PDF ###########################
dUTLD <- function(x, theta) {
  if (theta == 1) stop("theta cannot be exactly 1")
  pdf <- (log(theta) / (theta - 1)) * (theta^(x))
  return(pdf)
}


ll.UTLD <- function(delta)
{
  theta <- delta[1]
  
  if (theta <= 0 || theta == 1)
  {
    return(-Inf)
  }
  
  f <- dUTLD(x, theta)
  l <- log(f)
  logL <- sum(l)
  return(logL)
}

k <- maxLik(ll.UTLD,
            start  = c(0.5), # first data c(0.2, 0.5, 0.5),
            method = "BFGS")

s_e <- sqrt(abs(diag(solve(-k$hessian))))
k   <- k$estimate
k

p <- length(k)

theta <- k[1]

like <- ll.UTLD(k)
LL   <- -2 * like
AIC  <- LL + 2 * p
HQIC <- LL + 2 * p * log(log(n))
est_par1 <- cbind(k[1], s_e[1])
est_par1

cdf.UTLD <- function(x){
  cdf <- (theta^x -1) /(theta - 1)
  return(cdf)
}

UTLD_ks  <- ks.test(x, cdf.UTLD)
UTLD_cvm <- cvm.test(x, cdf.UTLD)
UTLD_ad  <- ad.test(x, cdf.UTLD)
UTLD_all <- cbind(
  LL, AIC, HQIC,
  UTLD_ks$statistic,  UTLD_ks$p.value,
  UTLD_cvm$statistic, UTLD_cvm$p.value,
  UTLD_ad$statistic,  UTLD_ad$p.value
)

UTLD_all

pUTLD <- function(x, theta){
  cdf <- (theta^x -1) /(theta - 1)
  return(cdf)
}

lwd=2
cols= c("black","red","magenta","orange", "purple", "springgreen",
        "blue", "sienna", "brown", "pink")
#cols=3:11
windows(5,5)
par(mar = c(4.5, 4.5, 1.1, 1.1)) #, omr=c(2,2,0,0))

#par(mfrow=c(3,3))

hist(x, breaks = 5, probability = T, col="#b3b3b3", xlim=c(0, 1),
     ylab="Density", main = "", xlab= "x")
curve(dUTLD(x, est_par1[1]), add = T, col=cols[1], lwd=lwd)
box()

##### PP plots for All models under comparison UPLD ======
windows(5,5)
#par(mfrow = c(3, 3))
par(mar = c(4.5, 4.5, 1.1, 1.1))
plot(ppoints(length(x)), sort(pUTLD(x, est_par1[1])),
     col = "blue",
     xlab = "Observed Probability", ylab = "Expected Probability", main = "")
abline(0,1)



library(AdequacyModel)

windows(8,4)
par(mfrow=c(1,2), mar=c(4.5,4.5,1.1,1))

# ---------------- TTT Plot ----------------
TTT(x,
    col  = "#00C98D",
    lwd  = 3,
    grid = TRUE,
    lty  = 2)

# title INSIDE plot region
usr <- par("usr")
text(x = usr[1] + 0.25*(usr[2]-usr[1]),
     y = usr[4] - 0.08*(usr[4]-usr[3]),
     labels = "TTT Plot",
     font = 1.5, cex = 1.0)

# ---------------- Box Plot ----------------
boxplot(x,
        horizontal = TRUE,
        col = "#4DBBD5",
        border = "#2C3E50",
        notch = TRUE,
        lwd = 2,
        boxwex = 0.5,
        staplewex = 0.7,
        outpch = 19,
        outcol = "#E64B35",
        xlab = "Data")

points(mean(x),1, pch=18, col="#D81B60", cex=1.6)
grid(nx=NA, ny=NULL, lty=3, col="gray80")

# title INSIDE plot region
usr <- par("usr")
text(x = usr[1] + 0.25*(usr[2]-usr[1]),
     y = usr[4] - 0.10*(usr[4]-usr[3]),
     labels = "Box Plot",
     font = 1.5, cex = 1.0)








