
library(goftest)
library(maxLik)

rm(list=ls())

x <-c(0.023,	0.032,	0.054,	0.069,	0.081,	0.094,
      0.105,	0.127,	0.148,	0.169,	0.188,	0.216,
      0.255,	0.277,	0.311,	0.361,	0.376,	0.395,
      0.432,	0.463,	0.481,	0.519,	0.529,	0.567,
      0.642,	0.674,	0.752,	0.823,	0.887,	0.926)

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
     font = 2, cex = 1.3)

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
     font = 2, cex = 1.3)




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

### Bounded Zeghdoudi distribution (BZD) ####
#### log likelihood function "BZD"
ll.BZD <- function(para) {
  epsilon <- para[1]
  if ( epsilon <= 0)
  {
    return (-Inf)
  }
  pdf <- (epsilon^3 * x^(epsilon - 1) * (log(1 / x) + 1) *
            log(1 / x)) / (epsilon + 2)
  logL <- sum(log(pdf))
  return(logL)
}


k<-maxLik(ll.BZD, start=c(1.0025), method="BFGS")
s_e<- sqrt(abs(diag(solve(-k$hessian))))
k<- k$estimate
k
p<- length(k) # no of parameters
epsilon<- k[1]
like<- ll.BZD(k)
LL<- -2*like
AIC<- LL + 2*p
HQIC <- LL + 2 * p * log(log(n))
est_par2<- cbind(k[1], s_e[1])

#============ CDF function ========

cdf<- function (x){
  cdf<-x^epsilon * ((epsilon^2 * log(x)^2 - epsilon *
       (epsilon + 2) * log(x)) / (epsilon + 2) + 1)
  return(cdf)
} 

BZD_ks<- ks.test(x, cdf)
BZD_cvm<- cvm.test(x, cdf)
BZD_ad<- ad.test(x, cdf)
BZD_all<- cbind(LL, AIC, HQIC, BZD_ks$statistic, BZD_ks$p.value,
                BZD_cvm$statistic, BZD_cvm$p.value, 
                BZD_ad$statistic, BZD_ad$p.value)
BZD_all

#### PDF of BZD
dBZD<- function(x, epsilon){
  pdf <- (epsilon^3 * x^(epsilon - 1) * (log(1 / x) + 1) *
            log(1 / x)) / (epsilon + 2)
  return(pdf)
}

#### CDF of BZD
pBZD<- function(x, epsilon){
  cdf <- x^epsilon * ((epsilon^2 * log(x)^2 - epsilon *
                         (epsilon + 2) * log(x)) / (epsilon + 2) + 1)
  return(cdf)
}


###Exponented Reduced Kies distribution ExRKD
ll.ExRKD <- function(para){
  beta <- para[1]
  delta <- para[2]
  u <- x/(1-x)
  logL <- sum(log(beta*delta)+(beta-1)*log(x)-u^beta-(1-x)^(beta+1)
              +(delta-1)*log(1-exp(-u^beta)))
  return(logL)
}

k<-maxLik(ll.ExRKD, start=c(0.5, 1.25), method="BFGS")
s_e<- sqrt(abs(diag(solve(-k$hessian))))
k<- k$estimate
k
p<- length(k) # no of parameters
beta <- k[1]
delta <- k[2]
like<- ll.ExRKD(k)
LL<- -2*like
AIC<- LL + 2*p
HQIC <- LL + 2 * p * log(log(n))
est_par3<- cbind(k[1], s_e[1], k[2], s_e[2])
#============ CDF function ========
cdf<- function (x){
  u <- x/(1-x)
  cdf <- (1-exp(-u^beta))^delta
  return(cdf)
}
ExRKD_ks<- ks.test(x, cdf)
ExRKD_cvm<- cvm.test(x, cdf)
ExRKD_ad<- ad.test(x, cdf)
ExRKD_all<- cbind(LL, AIC, HQIC, ExRKD_ks$statistic, ExRKD_ks$p.value,
                  ExRKD_cvm$statistic, ExRKD_cvm$p.value, 
                  ExRKD_ad$statistic, ExRKD_ad$p.value)
ExRKD_all

# PDF of Unit ExRKD 
dExRKD<-function(x, beta, delta){
  u<- (x/(1-x))^beta
  v <- exp(1-u)
  pdf<- beta*delta*x^(beta-1)*v*(1-v)^(delta-1)/((1-x)^(beta+1))
  return(pdf)
}
# CDF of Unit ExRKD 
pExRKD<-function(x, beta, delta){
  u <- x/(1-x)
  cdf <- (1-exp(-u^beta))^delta
  return(cdf)
}

## Kumarashwamy
#### log likelihood function "Log-LL"
ll.Ksm <- function(para) {
  alpha <- para[1]
  beta <- para[2]
  logL <- sum(log(alpha*beta)+(alpha-1)*log(x)+(beta-1)*log(1-x^alpha))
  return(logL)
}

k<-maxLik(ll.Ksm, start=c(0.25, 0.75), method="BFGS")
s_e<- sqrt(abs(diag(solve(-k$hessian))))
k<- k$estimate
k
p<- length(k) # no of parameters
alpha <- k[1]
beta <- k[2]
like <- ll.Ksm(k)
LL <- -2*like
AIC <- LL + 2*p
HQIC <- LL + 2 * p * log(log(n))
est_par4 <- cbind(k[1], s_e[1], k[2], s_e[2])
#============ CDF function ========
cdf<- function (x){
  z<- x^alpha
  cdf <- 1-(1-z)^beta
  return(cdf)
}

Ksm_ks<- ks.test(x, cdf)
Ksm_cvm<- cvm.test(x, cdf)
Ksm_ad<- ad.test(x, cdf)
Ksm_all<- cbind(LL, AIC, HQIC, Ksm_ks$statistic, Ksm_ks$p.value,
                Ksm_cvm$statistic, Ksm_cvm$p.value, 
                Ksm_ad$statistic, Ksm_ad$p.value)
Ksm_all

## PDF of Ksm
dKsm<- function(x, alpha, beta){
  u<- (1-x^alpha)^(beta-1)
  pdf<- alpha*beta*x^(alpha-1)*u
  return(pdf)
}

##cdf of Ksm
pKsm<- function(x, alpha, beta){
  z<- x^alpha
  cdf <- 1-(1-z)^beta
  return(cdf)
}

#### log likelihood function "beta"
ll.BETAD <- function(para){
  alpha <- para[1]
  beta <- para[2]
  logL <- sum((alpha - 1) * log(x) + (beta - 1) * 
                log(1 - x) - log(beta(alpha, beta)))
  return(logL)
}

k<-maxLik(ll.BETAD, start=c(0.25, 1.5), method="BFGS")
s_e<- sqrt(abs(diag(solve(-k$hessian))))
k<- k$estimate
k
p<- length(k) # no of parameters
alpha <- k[1]
beta <- k[2]
like<- ll.BETAD(k)
LL<- -2*like
AIC<- LL + 2*p
HQIC <- LL + 2 * p * log(log(n))
est_par5<- cbind(k[1], s_e[1], k[2], s_e[2])
#============ CDF function ========
cdf<- function(x) pbeta(x, alpha, beta)
BETAD_ks<- ks.test(x, cdf)
BETAD_cvm<- cvm.test(x, cdf)
BETAD_ad<- ad.test(x, cdf)
BETAD_all<- cbind(LL, AIC, HQIC, BETAD_ks$statistic, BETAD_ks$p.value,
                  BETAD_cvm$statistic, BETAD_cvm$p.value, 
                  BETAD_ad$statistic, BETAD_ad$p.value)
BETAD_all

#Unit Beta
dBETAD<-function(x, alpha, beta){
  pdf<- dbeta(x, alpha, beta)
  return(pdf)
}
### CDF Beta
pBETAD<-function(x, alpha, beta){
  cdf<- pbeta(x, alpha, beta)
  return(cdf)
}

#### log likelihood function "Unit Burr-III"
ll.UBIII <- function(para){
  beta <- para[1]
  lambda <- para[2]
  u<- (1/x)-1
  logL <- sum(log(beta*lambda*x^-2)+(beta-1)*log(u)-
                (lambda+1)*log(1+u^beta))
  return(logL)
}

k<-maxLik(ll.UBIII, start=c(0.25, 0.5), method="BFGS")
s_e<- sqrt(abs(diag(solve(-k$hessian))))
k<- k$estimate
k
p<- length(k) # no of parameters
beta <- k[1]
lambda <- k[2]
like<- ll.UBIII(k)
LL<- -2*like
AIC<- LL + 2*p
HQIC <- LL + 2 * p * log(log(n))
est_par6<- cbind(k[1], s_e[1], k[2], s_e[2])
#============ CDF function ========
cdf<- function (x){
  u<- (1/x)-1
  cdf <- (1+u^beta)^-lambda 
  return(cdf)
} 
UBIII_ks<- ks.test(x, cdf)
UBIII_cvm<- cvm.test(x, cdf)
UBIII_ad<- ad.test(x, cdf)
UBIII_all<- cbind(LL, AIC, HQIC, UBIII_ks$statistic, UBIII_ks$p.value,
                  UBIII_cvm$statistic, UBIII_cvm$p.value, 
                  UBIII_ad$statistic, UBIII_ad$p.value)
UBIII_all

### PDF of UBIII distribution

dUBIII<- function(x, beta, lambda){
  u<- (1/x)-1
  pdf<- beta*lambda*x^(-2)*u^(beta-1)*(1+u^beta)^(-lambda-1)
  return(pdf)
}

##cdf of UBIII
pUBIII<- function(x, beta, lambda){
  u<- (1/x)-1
  cdf <- (1+u^beta)^-lambda
}


### 6 ## Power Upper Truncated Weibull 
#### log likelihood function "PUTWD"
ll.PUTWD <- function(para){
  xi <- para[1]
  tau <- para[2]
  epsilon <- para[3]
  f_x<-  (xi * tau * epsilon * x^(epsilon - 1) * (1 - x^epsilon)^(tau - 1) * 
            exp(xi - xi * (1 - x^epsilon)^tau)) / (exp(xi) - 1)
  ll <- return(sum(log(f_x)))
}

k<-maxLik(ll.PUTWD, start = c(5.0, 0.5, 0.25), method="BFGS")
s_e <- sqrt(abs(diag(solve(-k$hessian))))
k <- k$estimate
k
p <- length(k) # no of parameters
xi <- k[1]
tau <- k[2]
epsilon <- k[3]
like<- ll.PUTWD(k)
LL<- -2*like
AIC<- LL + 2*p
#BIC <- LL + p * log(n)                  
#CAIC <- AIC+ (2*p*(p+1)/(n-p+1))
HQIC <- LL + 2 * p * log(log(n))
est_par7 <- cbind(k[1], s_e[1], k[2], s_e[2], k[3], s_e[3])

cdf<- function(x){
  cdf<- 1 - (1 - exp(-xi * (1 - x^epsilon)^tau)) / (1 - exp(-xi))
  return(cdf)
} 
PUTWD_ks <- ks.test(x, cdf)
PUTWD_cvm <- cvm.test(x, cdf)
PUTWD_ad <- ad.test(x, cdf)
PUTWD_all <- cbind(LL, AIC, HQIC, PUTWD_ks$statistic, PUTWD_ks$p.value,
                   PUTWD_cvm$statistic, PUTWD_cvm$p.value, 
                   PUTWD_ad$statistic, PUTWD_ad$p.value)
PUTWD_all

# Define the CDF function of PUTWD
pPUTWD <- function(x, xi, tau, epsilon) {
  cdf <- 1 - (1 - exp(-xi * (1 - x^epsilon)^tau)) / (1 - exp(-xi))
  return(cdf)
}

# Define the PDF function of PUTWD
dPUTWD <- function(x, xi, tau, epsilon) {
  pdf <- (xi * tau * epsilon * x^(epsilon - 1) * (1 - x^epsilon)^(tau - 1) * 
            exp(xi - xi * (1 - x^epsilon)^tau)) / (exp(xi) - 1)
  return(pdf)
}

#### log likelihood function "ExUTD"
ll.ExUTD <- function(para) {
  alpha <- para[1]
  beta <- para[2]
  if ( alpha<= 0 || beta < 1 )
  {
    return (-Inf)
  }  
  logL <- sum(log(alpha)-(2*alpha+1)*log(x)+ 
                beta-beta*x^(-alpha)+log(beta-x^(alpha)))
  return(logL)
}


k<-maxLik(ll.ExUTD, start=c(0.5, 2.25), method="BFGS") # beta>1
s_e<- sqrt(abs(diag(solve(-k$hessian))))
k<- k$estimate
k
p<- length(k) # no of parameters
alpha<- k[1]
beta<- k[2]
like<- ll.ExUTD(k)
LL<- -2*like
AIC<- LL + 2*p
HQIC <- LL + 2 * p * log(log(n))
est_par8<- cbind(k[1], s_e[1], k[2], s_e[2])
#============ CDF function ========
cdf<- function (x){
  cdf <- x^(-alpha)*exp(beta*(1-x^(-alpha)))  
} 
ExUTD_ks<- ks.test(x, cdf)
ExUTD_cvm<- cvm.test(x, cdf)
ExUTD_ad<- ad.test(x, cdf)
ExUTD_all<- cbind(LL, AIC, HQIC, ExUTD_ks$statistic, ExUTD_ks$p.value,
                  ExUTD_cvm$statistic, ExUTD_cvm$p.value, 
                  ExUTD_ad$statistic, ExUTD_ad$p.value)
ExUTD_all

# PDF of ExUTD
dExUTD <- function (x, alpha, beta)
{
  u<- x^(-alpha)
  pdf <- alpha*(x^(-2*alpha-1))*exp(beta-beta*u)*(beta-x^(alpha))
  return(pdf)
}

# CDF of ExUTD
pExUTD <- function (x, alpha, beta)
{
  cdf <- x^(-alpha)*exp(beta*(1-x^(-alpha)))
  return(cdf)
}

#### log likelihood function "ExTL"
ll.ExTL <- function(para) {
  alpha <- para[1]
  beta <- para[2]
  lambda <- para[3]
  u <- x/beta
  v<- 1-(u^alpha)*(2-u)^alpha
  logL <- sum(log(2*alpha*lambda/beta)+log(1-u)+(alpha-1)*log(u*(2-u))+
                (lambda-1)*log(v))
  return(logL)
}


k<-maxLik(ll.ExTL, start=c(0.25, 1.5, 1.25), method="BFGS")
s_e<- sqrt(abs(diag(solve(-k$hessian))))
k<- k$estimate
k
p<- length(k) # no of parameters
alpha <- k[1]
beta <- k[2]
lambda <- k[3]
like<- ll.ExTL(k)
LL<- -2*like
AIC<- LL + 2*p
HQIC <- LL + 2 * p * log(log(n))
est_par9<- cbind(k[1], s_e[1], k[2], s_e[2], k[3], s_e[3])
#============ CDF function ========
cdf<- function (x){
  u <- x/beta
  v<- 1-(u^alpha)*(2-u)^alpha
  cdf<- 1- v^lambda
  return(cdf)
} 
ExTL_ks<- ks.test(x, cdf)
ExTL_cvm<- cvm.test(x, cdf)
ExTL_ad<- ad.test(x, cdf)
ExTL_all<- cbind(LL, AIC, HQIC, ExTL_ks$statistic, ExTL_ks$p.value,
                 ExTL_cvm$statistic, ExTL_cvm$p.value, 
                 ExTL_ad$statistic, ExTL_ad$p.value)
ExTL_all

## PDF of ExTL *********************************
dExTL<- function(x, alpha, beta, lambda){
  u<- x/beta
  pdf<- (2*alpha*lambda/beta)*(1-u)*(u*(2-u))^(alpha-1)*
    (1-u^(alpha)*(2-u)^(alpha))^(lambda-1)
  return(pdf)
}
#### CDF of ExTL  ####
pExTL<- function(x, alpha, beta, lambda){
  u <- x/beta
  v<- 1-(u^alpha)*(2-u)^alpha
  cdf<- 1- v^lambda
  return(cdf)
}


#-----------------------------
# Unit Topp-Leone Distribution
#-----------------------------

# PDF
dTLD <- function(x, alpha){
  ifelse(x>0 & x<1,
         2*alpha*(1-x)*(2*x - x^2)^(alpha-1),
         0)
}

# CDF
pTLD <- function(x, alpha){
  ifelse(x>0 & x<1,
         (2*x - x^2)^alpha,
         ifelse(x>=1,1,0))
}

ll.TLD <- function(delta)
{
  alpha <- delta[1]
  f <- dTLD(x, alpha)
  l <- log(f)
  logL <- sum(l)
  return(logL)
}

k <- maxLik(ll.TLD,
            start  = c(0.05), # first data c(0.2, 0.5, 0.5),
            method = "BFGS")

s_e <- sqrt(abs(diag(solve(-k$hessian))))
k   <- k$estimate
k

p <- length(k)

alpha <- k[1]

like <- ll.TLD(k)
LL   <- -2 * like
AIC  <- LL + 2 * p
HQIC <- LL + 2 * p * log(log(n))
est_par10 <- cbind(k[1], s_e[1])
est_par10

cdf.TLD <- function(x){
  cdf <-  (2*x - x^2)^alpha
  return(cdf)
}

TLD_ks  <- ks.test(x, cdf.TLD)
TLD_cvm <- cvm.test(x, cdf.TLD)
TLD_ad  <- ad.test(x, cdf.TLD)
TLD_all <- cbind(
  LL, AIC, HQIC,
  TLD_ks$statistic,  TLD_ks$p.value,
  TLD_cvm$statistic, TLD_cvm$p.value,
  TLD_ad$statistic,  TLD_ad$p.value
)

TLD_all


#-----------------------------
# Unit Lindley Distribution
#-----------------------------

# PDF
dULD <- function(x, theta){
  pdf <- (theta^2/(1+theta))*(1/(1-x)^3)*
    exp(-theta*x/(1-x))
  return(pdf)
}

# CDF
pULD <- function(x, theta){
  cdf <- 1 - (1 - theta*x/((1+theta)*(x-1)))*
    exp(-theta*x/(1-x))
  return(cdf)
}


ll.ULD <- function(para) {
  theta <- para[1]
  if (theta <= 0)
  {
    return (-Inf)
  }
  pdf <-  (theta^2/(1+theta))*(1/(1-x)^3)*
    exp(-theta*x/(1-x))
  logL <- sum(log(pdf))
  return(logL)
}


k<-maxLik(ll.ULD, start=c(0.05), method="BFGS")
s_e<- sqrt(abs(diag(solve(-k$hessian))))
k<- k$estimate
k
p<- length(k) # no of parameters
theta<- k[1]
like<- ll.ULD(k)
LL<- -2*like
AIC<- LL + 2*p
HQIC <- LL + 2 * p * log(log(n))
est_par11<- cbind(k[1], s_e[1])

#============ CDF function ========

cdf<- function (x){
  cdf1<-  1 - (1 - theta*x/((1+theta)*(x-1)))*
    exp(-theta*x/(1-x))
  return(cdf1)
} 

ULD_ks<- ks.test(x, cdf)
ULD_cvm<- cvm.test(x, cdf)
ULD_ad<- ad.test(x, cdf)
ULD_all<- cbind(LL, AIC, HQIC, ULD_ks$statistic, ULD_ks$p.value,
                ULD_cvm$statistic, ULD_cvm$p.value, 
                ULD_ad$statistic, ULD_ad$p.value)
ULD_all


#-----------------------------
# Unit Teissier Distribution
#-----------------------------

# PDF
dUTeD <- function(x, theta){
  pdf <- theta*(1/x^theta -1)*x^(-(theta+1))*
    exp((-x^(-theta))+1)
  return(pdf)
}

# CDF
pUTeD <- function(x, theta){
  cdf <- x^(-theta)*exp((-x^(-theta))+1)
  return(cdf)
}


ll.UTeD <- function(para) {
  theta <- para[1]
  if (theta <= 0)
  {
    return (-Inf)
  }
  pdf <- theta*(1/x^theta -1)*x^(-(theta+1))*
    exp((-x^(-theta))+1)
  logL <- sum(log(pdf))
  return(logL)
}


k<-maxLik(ll.UTeD, start=c(0.05), method="BFGS")
s_e<- sqrt(abs(diag(solve(-k$hessian))))
k<- k$estimate
k
p<- length(k) # no of parameters
theta<- k[1]
like<- ll.UTeD(k)
LL<- -2*like
AIC<- LL + 2*p
HQIC <- LL + 2 * p * log(log(n))
est_par12<- cbind(k[1], s_e[1])

#============ CDF function ========

cdf<- function (x){
  cdf1<- x^(-theta)*exp((-x^(-theta))+1)
  return(cdf1)
} 

UTeD_ks<- ks.test(x, cdf)
UTeD_cvm<- cvm.test(x, cdf)
UTeD_ad<- ad.test(x, cdf)
UTeD_all<- cbind(LL, AIC, HQIC, UTeD_ks$statistic, UTeD_ks$p.value,
                 UTeD_cvm$statistic, UTeD_cvm$p.value, 
                 UTeD_ad$statistic, UTeD_ad$p.value)
UTeD_all


#-----------------------------
# Unit Rayleigh Distribution
#-----------------------------

# PDF
dURayD <- function(x, alpha){
  pdf<- (-2*alpha/x) *log(x)*exp(- alpha* (log(x))^2)
  return(pdf)
}

# CDF
pURayD <- function(x, alpha){
  cdf<-exp(- alpha* (log(x))^2)
  return(cdf)
}

ll.URayD <- function(delta)
{
  alpha <- delta[1]
  f <- dURayD(x, alpha)
  l <- log(f)
  logL <- sum(l)
  return(logL)
}

k <- maxLik(ll.URayD,
            start  = c(0.5), # first data c(0.2, 0.5, 0.5),
            method = "BFGS")

s_e <- sqrt(abs(diag(solve(-k$hessian))))
k   <- k$estimate
k

p <- length(k)

alpha <- k[1]

like <- ll.URayD(k)
LL   <- -2 * like
AIC  <- LL + 2 * p
HQIC <- LL + 2 * p * log(log(n))
est_par13 <- cbind(k[1], s_e[1])
est_par13

cdf.URayD <- function(x){
  cdf <- exp(- alpha* (log(x))^2)
  return(cdf)
}

URayD_ks  <- ks.test(x, cdf.URayD)
URayD_cvm <- cvm.test(x, cdf.URayD)
URayD_ad  <- ad.test(x, cdf.URayD)
URayD_all <- cbind(
  LL, AIC, HQIC,
  URayD_ks$statistic,  URayD_ks$p.value,
  URayD_cvm$statistic, URayD_cvm$p.value,
  URayD_ad$statistic,  URayD_ad$p.value
)

URayD_all




#%%%%%%%%%%%% Data and table Creation %%%%%%%%%%%%%%%
library(dplyr)

est_par_se <- bind_rows(as.data.frame(est_par1), as.data.frame(est_par2),
                         as.data.frame(est_par4),
                        as.data.frame(est_par5), as.data.frame(est_par6),
                         as.data.frame(est_par10),
                        as.data.frame(est_par11), as.data.frame(est_par12)
                        , as.data.frame(est_par13))
rownames(est_par_se) <- c("UPLD", "BZD", "Ksm", 
                          "BETAD", "UBIII", "TLD", "ULD", "UTeD", "URayD")
colnames(est_par_se) <- c("parameter", "SE", "parameter", "SE")
est_par_se<- round(est_par_se, 4)
print(est_par_se)

#$$$$$$$$$$$$$$$$$$$$ Table for AIC....%%%%%%%%%%%%%%%%%%%%

stat.aic<- bind_rows(as.data.frame(UTLD_all),
                     as.data.frame(BZD_all), 
                     as.data.frame(Ksm_all), as.data.frame(BETAD_all),
                     as.data.frame(UBIII_all),
                     as.data.frame(TLD_all), as.data.frame(ULD_all),
                     as.data.frame(UTeD_all),  as.data.frame(URayD_all))
rownames(stat.aic)<- c("UPLD", "BZD", "Ksm", 
                       "BETAD", "UBIII", "TLD", "ULD", "UTeD", "URayD")
colnames(stat.aic)<- c("-2logL", "AIC", "HQIC", "KS", "p(KS)",
                       "CVM", "p(CVM)", "AD", "p(AD)")

stat.aic<- round(stat.aic, 4)
stat.aic



write.csv(est_par_se, "MLEs_se_UPLD_set1.csv", row.names=T)

write.csv(stat.aic, "AIC_KS_table_UPLD_set1.csv", row.names=T)


lwd=2
cols= c("black","red","magenta","orange", "purple", "springgreen",
        "blue", "sienna", "brown", "pink")
#cols=3:11
windows(7,7)
par(mar = c(4.5, 4.5, 1.1, 1.1), omr=c(2,2,0,0))

par(mfrow=c(3,3))

hist(x, breaks = 10, probability = T, col="#b3b3b3", xlim=c(0, 1),
      ylab="Density", main = "UPLD", xlab= "x")
curve(dUTLD(x, est_par1[1]), add = T, col=cols[1], lwd=lwd)  #ylim=c(0, 7.5),

#hist(x, breaks = 10, probability = T, col="#b3b3b3", xlim=c(0, 0.6),
#     ylim=c(0, 7.5), ylab="Density", main = "UTMAD", xlab= "x")
#curve(dUTMAD(x, est_par2[1], est_par2[3]), add=T, col=cols[2], lwd=lwd)


hist(x, breaks = 10, probability = T, col="#b3b3b3", xlim=c(0, 1),
      ylab="Density", main = "BZD", xlab= "x")
curve(dBZD(x, est_par2[1]), add=T, col=cols[3],lwd=lwd)

#hist(x, breaks = 10, probability = T, col="#b3b3b3", xlim=c(0, 1),
#      ylab="Density", main = "ExRKD", xlab= "x")
#curve(dExRKD(x, est_par3[1], est_par3[3]), add=T, col=cols[4],lwd=lwd)

hist(x, breaks = 10, probability = T, col="#b3b3b3", xlim=c(0, 1),
      ylab="Density", main = "KsmD", xlab= "x")
curve(dKsm(x, est_par4[1], est_par4[3]), add=T, col=cols[5],lwd=lwd)

hist(x, breaks = 10, probability = T, col=  "#b3b3b3", xlim=c(0, 1),
      ylab="Density", main = "BetaD", xlab= "x")
curve(dBETAD(x, est_par5[1], est_par5[3]), add=T, col=cols[6],lwd=lwd)

hist(x, breaks = 10, probability = T, col=  "#b3b3b3", xlim=c(0, 1),
     ylab= "Density", main = "UBIIID", xlab= "x")
curve(dUBIII(x, est_par6[1], est_par6[3]), add=T, col=cols[7],lwd=lwd)

#hist(x, breaks = 10, probability = T, col=  "#b3b3b3", xlim=c(0, 1),
#      ylab= "Density", main = "PUTWD", xlab= "x")
#curve(dPUTWD(x, est_par7[1], est_par7[3], est_par7[5]), add=T, col=cols[8],lwd=lwd)

#hist(x, breaks = 10, probability = T, col=  "#b3b3b3", xlim=c(0, 1),
#      ylab= "Density", main = "ExUTD", xlab= "x")
#curve(dExUTD(x, est_par8[1], est_par8[3]), add=T, col=cols[9],lwd=lwd)

#hist(x, breaks = 10, probability = T, col=  "#b3b3b3", xlim=c(0, 1),
#      ylab= "Density", main = "ExTLD", xlab= "x")
#curve(dExTL(x, est_par9[1], est_par9[3], est_par9[5]), add=T, col=cols[10],lwd=lwd)

hist(x, breaks = 10, probability = T, col="#b3b3b3", xlim=c(0, 1),
     ylab="Density", main = "TLD", xlab= "x")
curve(dTLD(x, est_par10[1]), add=T, col=cols[2],lwd=lwd)

hist(x, breaks = 10, probability = T, col="#b3b3b3", xlim=c(0, 1),
      ylab="Density", main = "ULD", xlab= "x")
curve(dULD(x, est_par11[1]), add=T, col=cols[4],lwd=lwd)

hist(x, breaks = 10, probability = T, col="#b3b3b3", xlim=c(0, 1),
     ylim=c(0, 4),    ylab="Density", main = "UTeD", xlab= "x")
curve(dUTeD(x, est_par12[1]), add=T, col=cols[5],lwd=lwd)

hist(x, breaks = 10, probability = T, col=  "#b3b3b3", xlim=c(0, 1),
     ylim=c(0, 3),    ylab= "Density", main = "URayD", xlab= "x")
curve(dURayD(x, est_par13[1]), add=T, col=cols[10],lwd=lwd)


##### PP plots for All models under comparison UPLD ======
windows(7,7)
par(mfrow = c(3, 3))
par(mar = c(4.5, 4.5, 2.1, 1.1))
plot(ppoints(length(x)), sort(pUTLD(x, est_par1[1])),
     col = "blue",
     xlab = "", ylab = "Expected Probability", main = "UPLD")
abline(0,1)

plot(ppoints(length(x)), sort(pBZD(x, est_par2[1])), col = "red",
     xlab = "", ylab = "", main = "BZD")
abline(0,1)

plot(ppoints(length(x)), sort(pKsm(x, est_par4[1], est_par4[3])),
     col = "purple", xlab = "", ylab = "", main = "KsmD")
abline(0,1)

plot(ppoints(length(x)), sort(pBETAD(x, est_par5[1], est_par5[3])), col = "gray",
     xlab = "", ylab = "Expected Probability", main = "BETAD")
abline(0,1)

plot(ppoints(length(x)), sort(pUBIII(x, est_par6[1], est_par6[3])),
     col = "orange", xlab = "", ylab = "", main = "UBIIID")
abline(0,1)

plot(ppoints(length(x)), sort(pTLD(x, est_par10[1])), col = "green",
     xlab = "", ylab = "", main = "TLD")
abline(0,1)

plot(ppoints(length(x)), sort(pULD(x, est_par11[1])),
     col = "magenta", xlab = "Observed Probability", 
     ylab = "Expected Probability", main = "ULD")
abline(0,1)

plot(ppoints(length(x)), sort(pUTeD(x, est_par12[1])), col = "black",
     xlab = "Observed Probability", ylab = "", main = "UTeD")
abline(0,1)

plot(ppoints(length(x)), sort(pURayD(x, est_par13[1])),
     col = "brown", xlab = "Observed Probability", ylab = "", main = "URayD")
abline(0,1)




# The Food Insecurity Multidimensional Index(FIMI) Data good for all
x<- c(0.0666, 0.4524, 0.0878,0.5505,0.1236,0.7038,0.1964,0.7569,0.4212,
      0.7754,0.4378,0.9608,0.442,0.9795)

x<-c(0.0009, 0.0040, 0.0142, 0.0221, 0.0261, 0.0418, 0.0473, 0.0834, 
     0.1091, 0.1252, 0.1404, 0.1498, 0.1750, 0.2031, 0.2099, 0.2168,
     0.2918, 0.3465, 0.4035, 0.6143)

#rain_prop
x <- c(
  0.00,0.03,0.10,0.25,0.40,0.63,0.83,0.93,1.00,0.97,
  0.87,0.70,0.52,0.35)



url <- "https://data.weather.gov.hk/weatherAPI/cis/csvfile/PLC/ALL/daily_PLC_RF_ALL.csv"

rain <- read.csv(
  url,
  skip = 2,          # skip metadata rows at top
  fill = TRUE,       # fill uneven rows
  header = TRUE,
  stringsAsFactors = FALSE
)

names(rain) <- c("Year","Month","Day","Rainfall_mm","Data_Completeness")
head(rain, 50)

rain$Date <- as.Date(with(rain, paste(Year, Month, Day, sep = "-")))
head(rain)

rain <- subset(rain, Data_Completeness == "C")

rain$Rainfall_mm <- as.numeric(rain$Rainfall_mm)

rain2026 <- subset(rain, Year == 2026)

x <- c(16, 54.5, 83.5, 64.5, 30.4, 34.5, 0.0, 0.0, 0.0, 0.0, 15.0, 2.5, 0.0, 
       78.0, 2.5, 0.0, 36.5, 113.5, 87.5, 3.5, 0.0, 0.0, 0.5, 0.0, 2.0, 0.0,
       0.0, 0.7, 0.0,
       0.5, 4.5)
x<-x/max(x)
hist(x)

rain <- na.omit(rain)

# convert to [0,1]
rain_prop <- rain2026$Rainfall_mm

# remove NA
rain_prop <- rain_prop/ max(rain_prop)
x <- rain_prop
summary(rain_prop) 




# Install the package if you haven't already
install.packages("betareg")

# Load the dataset
library(betareg)
data("FoodExpenditure")

# Check the range of the proportion variable
range(FoodExpenditure$food)
# Note: This may be within (0,1) and not include the boundaries


library(betareg)
data("ReadingSkills")
range(ReadingSkills$accuracy)   # Check: includes 0 and 1?
# This dataset is often used for beta regression with exact 0/1

x<- ReadingSkills$accuracy1


# Download and load the AlcoholUse dataset directly from CRAN
load(url("https://github.com/cran/zoib/raw/master/data/AlcoholUse.rda"))

# Now the dataset is available as an object named "AlcoholUse"
head(AlcoholUse)

# Verify it includes 0 and 1
range(AlcoholUse$Percentage)

x<- AlcoholUse$Percentage

install.packages("geepack")
library(geepack)
data("dietox")
# Create a proportion variable, e.g., weight gain relative to final weight
# This may require transformation, but the dataset itself is real.


