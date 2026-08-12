

library(goftest)
library(maxLik)
#This is the data for pass rate of students in SLS in various subjects from
#Far Western Region in 2071BS

x <- c(52.04,79.81,45.88,46.70,0.00,70.52,88.47,87.18,50.00,77.98,96.88,92.82,
       78.31,63.64,100.00,50.00,98.71,94.78,88.67,78.05,86.77,93.73,100.00,
       88.55)
x<- x/100

n <- length(x)

library(AdequacyModel)

windows(8,4)
par(mfrow=c(1,2), mar=c(4.5,4.5,1.5,1))

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

points(mean(x),1, pch=18, col="#D81B60", cex=3)
grid(nx=NA, ny=NULL, lty=3, col="gray80")

# title INSIDE plot region
usr <- par("usr")
text(x = usr[1] + 0.25*(usr[2]-usr[1]),
     y = usr[4] - 0.10*(usr[4]-usr[3]),
     labels = "Box Plot",
     font = 1.5, cex = 1.0)




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

hist(x, breaks = 10, probability = T, col="#b3b3b3", xlim=c(0, 1),
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



Data 2
     Par       SE
 39.28543 23.72657

   LL       AIC      HQIC      KS         P-value     CVM         p-value AD p-value                        
-20.30449 -18.30449 -17.99196 0.1052464 0.9530817 0.04715371 0.8976081 Inf 2.5e-05


Data 3
   Par         SE
0.01106352 0.008658316   
   
LL         AIC      HQIC        KS        P-value    CVM      P-value   AD  P-value                            
-48.7192 -46.7192 -46.06973 0.1269233 0.4926062 0.1581136 0.3668332 Inf 1.395349e-05   
   
   