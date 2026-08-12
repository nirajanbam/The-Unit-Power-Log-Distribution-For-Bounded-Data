

library(goftest)
library(maxLik)
#This is the data for pass rate of students in SLS in Palpa

y <- c(
  0.2955, 0.3548, 0.3333, 0.0769, 0.3333, 0.3077, 0.3889, 0.3529, 0.2927, 0.3070, 0.4248,
  0.2281, 0.5849, 0.4583, 0.2182, 0.3380, 0.5714, 0.4839, 0.6250, 0.1778, 0.7917, 0.2791,
  0.4667, 0.1528, 0.1897, 0.2308, 0.2453, 0.1429, 0.3125, 0.3061, 0.4194, 0.2250, 0.4943,
  0.3111, 0.3023, 0.2500, 0.2958, 0.4250, 0.4074, 0.8333, 0.2653, 0.1296, 0.3056, 0.3725,
  0.2703, 0.3265, 0.6212, 0.2353, 0.0968, 0.4444, 0.3333, 0.4286, 0.9184, 0.2941, 0.0606,
  0.1556, 0.1176, 0.4333, 0.6471, 0.5000, 0.3056, 0.3462, 0.1750, 0.4000, 0.2069, 0.1250,
  0.1250, 0.5000, 0.8889, 0.2273, 0.2333, 0.4595, 0.2381, 0.2222, 0.0000, 0.4043, 0.5385,
  1.0000, 0.8235, 0.8261, 0.4359, 0.8519, 0.2000, 0.8182, 0.4500, 0.3684, 0.0250, 0.0833,
  0.1000, 0.9130, 0.2000, 0.2000, 0.4167, 0.5000, 0.4375, 0.1923, 0.5000, 0.9756, 0.1579,
  0.2000, 0.2955, 0.4048, 0.5556, 0.6667, 0.6111, 1.0000, 0.5000, 0.2500, 0.1000, 0.9000,
  0.4400, 0.3750, 0.5405, 0.3500, 0.5714, 0.3333, 0.4615
)  ## palpa
hist(y)
x<- y[75:105]

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
     col = "blue", pch=18,
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
   


Data2 Palpa
[1,] 0.5919873 0.3705721

       LL      AIC     HQIC      KS          P-value     CVM       P       AD           
D -0.7052084 1.294792 1.762236 0.09034147 0.9619939 0.05330129 0.8595356 Inf

D 1.935484e-05





  