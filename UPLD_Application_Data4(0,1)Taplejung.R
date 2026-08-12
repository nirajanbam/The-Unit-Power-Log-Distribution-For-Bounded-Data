
library(goftest)
library(maxLik)
# Data set 4 This is another data of proportion of students passed from each school in SLC
# from Taplejung in 2072 BS. 
x <- c(
  0.1186, 0.3448, 0.2041, 0.2917, 0.1429, 0.0959, 0.0690, 0.1579, 0.1731, 0.1538,
  0.1220, 0.1163, 0.2245, 0.0408, 0.0769, 0.2273, 0.4286, 0.0909, 0.2273, 1.0000,
  0.1471, 0.0400, 0.0667, 0.0930, 0.2667, 0.1538, 0.2895, 0.0811, 0.1622, 1.0000,
  0.0000, 0.8333, 0.1667, 0.0408, 0.1667, 0.7000, 0.0417, 0.2381, 0.1053, 0.1667,
  0.0000, 0.0000, 0.0000, 0.5714, 0.1429, 0.0116, 0.1250, 0.1500, 0.0645, 0.0909,
  0.0000, 0.1333, 0.0000, 0.0000, 0.0000, 0.1111, 0.0000, 0.0714, 0.0000, 0.0323,
  0.0000, 0.2667
)

x<-  x[1:32]        # x[1:32] # p-value 0.5173 # x[1:43]
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






