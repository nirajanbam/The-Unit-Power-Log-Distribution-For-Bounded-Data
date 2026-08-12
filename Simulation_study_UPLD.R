###  Simulation Experimennt ============ 
library(maxLik)
rm(list = ls())
set.seed(2083)

############################################################
# 1. Random Generator (Correct)
############################################################
rUPLD <- function(n, theta){
  u <- runif(n)
  log(1 + u*(theta-1))/log(theta)
}

############################################################
# 2. TRUE CDF (Corrected)
############################################################
pUPLD <- function(x, theta){
  (theta^x - 1)/(theta - 1)
}

############################################################
# 3. Method of Moments Estimator (Correct)
############################################################
MoM_UPLD <- function(x){
  xbar <- mean(x)
  
  f <- function(theta){
    xbar*(theta-1)*log(theta) - theta*log(theta) + theta - 1
  }
  
  uniroot(f, interval=c(1.01,50))$root
}

############################################################
# 4. Least Squares Estimator (Correct CDF used)
############################################################
LSE_UPLD <- function(x){
  n <- length(x)
  x <- sort(x)
  
  obj <- function(theta){
    Fth <- pUPLD(x, theta)
    Fi  <- (1:n)/(n+1)
    sum((Fth - Fi)^2)
  }
  
  optim(2, obj, method="L-BFGS-B", lower=1.01, upper=50)$par
}

############################################################
# 5. Simulation setup
############################################################
sample_sizes <- c(45,100,145,200,400,545,745,1045)
k <- 1000
theta_true <- 4.0 # set6=4.0   #set5=3.0  #set4=2.5   #set3=2.0  #set2=1.5 # set1=1.25
z <- qnorm(0.975)

Results <- data.frame()

############################################################
# 6. Monte Carlo Simulation
############################################################
for(n in sample_sizes){
  
  MLE <- MoM <- LSE <- SE_MLE <- numeric(k)
  
  for(i in 1:k){
    
    x <- rUPLD(n, theta_true)
    
    ################ MLE (ANALYTIC LOG-LIKELIHOOD) ################
    logLike <- function(theta){
      if(theta <= 0 || theta == 1) return(-Inf)
      n*log(log(theta)/(theta-1)) + log(theta)*sum(x)
    }
    
    fit <- tryCatch(
      maxLik(logLike,start=c(theta=2),method="BFGS"),
      error=function(e) NULL)
    
    if(!is.null(fit)){
      MLE[i] <- fit$estimate
      SE_MLE[i] <- sqrt(vcov(fit))
    } else {
      MLE[i] <- NA; SE_MLE[i] <- NA
    }
    
    ################ MoM ################
    MoM[i] <- tryCatch(MoM_UPLD(x), error=function(e) NA)
    
    ################ LSE ################
    LSE[i] <- tryCatch(LSE_UPLD(x), error=function(e) NA)
  }
  
  ## remove NA runs
  MLE <- na.omit(MLE)
  MoM <- na.omit(MoM)
  LSE <- na.omit(LSE)
  SE_MLE <- na.omit(SE_MLE)
  
  ###################################################
  # Monte-Carlo SE for MoM and LSE
  ###################################################
  SE_MoM <- sd(MoM)
  SE_LSE <- sd(LSE)
  
  ###################################################
  # Coverage probabilities (95% Wald-type)
  ###################################################
  CP_MLE <- mean(theta_true >= (MLE - z*SE_MLE) &
                   theta_true <= (MLE + z*SE_MLE))
  
  CP_MoM <- mean(theta_true >= (MoM - z*SE_MoM) &
                   theta_true <= (MoM + z*SE_MoM))
  
  CP_LSE <- mean(theta_true >= (LSE - z*SE_LSE) &
                   theta_true <= (LSE + z*SE_LSE))
  
  ###################################################
  # Performance measures
  ###################################################
  temp <- data.frame(
    SampleSize = n,
    
    Bias_MLE = mean(MLE) - theta_true,
    MSE_MLE  = mean((MLE-theta_true)^2),
    CP_MLE   = CP_MLE,
    
    Bias_MoM = mean(MoM) - theta_true,
    MSE_MoM  = mean((MoM-theta_true)^2),
    CP_MoM   = CP_MoM,
    
    Bias_LSE = mean(LSE) - theta_true,
    MSE_LSE  = mean((LSE-theta_true)^2),
    CP_LSE   = CP_LSE
  )
  
  Results <- rbind(Results,temp)
}

############################################################
# 7. Final Output
############################################################
Results <- round(Results,4)
print(Results)



write.csv(Results,
          "Simulation_UPLD_MLE_MoM_LSE_set6.csv",
         row.names = FALSE)

