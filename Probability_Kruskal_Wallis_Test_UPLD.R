library(maxLik)

#-----------------------------
# UPLD PDF & CDF
#-----------------------------
dUPLD <- function(x, theta){
  (log(theta)/(theta-1)) * theta^x
}

pUPLD <- function(x, theta){
  (theta^x - 1)/(theta - 1)
}

#-----------------------------
# Log-likelihood (must receive x!)
#-----------------------------
logLik_UPLD <- function(param, x)
{
  theta <- param[1]
  
  # parameter restriction
  if(theta <= 0 || abs(theta-1) < 1e-8)
    return(-Inf)
  
  f <- dUPLD(x, theta)
  
  # avoid log(0)
  if(any(f <= 0)) return(-Inf)
  
  sum(log(f))
}

#-----------------------------
# MLE function (FINAL)
#-----------------------------
fit_theta <- function(x)
{
  fit <- maxLik(logLik_UPLD,
                start = c(theta = 2),
                method = "BFGS",
                x = x)        # <-- pass data here
  
  coef(fit)[1]
}

data2 <- c(0.0000, 0.4043, 0.5385, 1.0000, 0.8235, 0.8261, 0.4359, 0.8519, 0.2000,
     0.8182, 0.4500, 0.3684, 0.0250, 0.0833, 0.1000, 0.9130, 0.2000, 0.2000,
     0.4167, 0.5000, 0.4375, 0.1923, 0.5000, 0.9756, 0.1579, 0.2000, 0.2955,
     0.4048, 0.5556, 0.6667, 0.6111)

data3 <- c(0.5935,0.6623,0.2917,0.7705,0.5000,0.4839,0.7536,0.9048,0.8667,0.9189,
           0.9118,0.7660,0.5111,0.2162,0.8485,0.8471,0.9091,0.2308,0.7209,0.9000,
           0.8000,0.4828,0.3333,0.8605,0.1250,1.0000,0.0000,0.0000,0.9870)

data4 <- c(0.1186,0.3448,0.2041,0.2917,0.1429,0.0959,0.0690,0.1579,0.1731,0.1538,
           0.1220,0.1163,0.2245,0.0408,0.0769,0.2273,0.4286,0.0909,0.2273,1.0000,
           0.1471,0.0400,0.0667,0.0930,0.2667,0.1538,0.2895,0.0811,0.1622,1.0000,
           0.0000,0.8333)


theta2 <- fit_theta(data2)
theta3 <- fit_theta(data3)
theta4 <- fit_theta(data4)

theta2; theta3; theta4

breaks <- c(0,0.01, seq(0.20,0.90,by=0.10),1.00)

interval_prob <- function(theta){
  probs <- numeric(length(breaks)-1)
  for(i in 1:(length(breaks)-1)){
    probs[i] <- pUPLD(breaks[i+1],theta) - pUPLD(breaks[i],theta)
  }
  probs
}

prob2 <- interval_prob(theta2)
prob3 <- interval_prob(theta3)
prob4 <- interval_prob(theta4)

interval_labels <- c("0-0.01","0.01-0.20","0.20-0.30","0.30-0.40","0.40-0.50",
                     "0.50-0.60","0.60-0.70","0.70-0.80","0.80-0.90","0.90-1.00")

result <- data.frame(
  Interval = interval_labels,
  DataII = round(prob2,4),
  DataIII = round(prob3,4),
  DataIV = round(prob4,4)
)

result


write.csv(result, "Computed_Probability.csv", row.names = FALSE)


# Kruskal-Wallis test
# combine data

data2 <- c(0.0000, 0.4043, 0.5385, 1.0000, 0.8235, 0.8261, 0.4359, 0.8519, 0.2000,
           0.8182, 0.4500, 0.3684, 0.0250, 0.0833, 0.1000, 0.9130, 0.2000, 0.2000,
           0.4167, 0.5000, 0.4375, 0.1923, 0.5000, 0.9756, 0.1579, 0.2000, 0.2955,
           0.4048, 0.5556, 0.6667, 0.6111)

data3 <- c(0.5935,0.6623,0.2917,0.7705,0.5000,0.4839,0.7536,0.9048,0.8667,0.9189,
           0.9118,0.7660,0.5111,0.2162,0.8485,0.8471,0.9091,0.2308,0.7209,0.9000,
           0.8000,0.4828,0.3333,0.8605,0.1250,1.0000,0.0000,0.0000,0.9870)

data4 <- c(0.1186,0.3448,0.2041,0.2917,0.1429,0.0959,0.0690,0.1579,0.1731,0.1538,
           0.1220,0.1163,0.2245,0.0408,0.0769,0.2273,0.4286,0.0909,0.2273,1.0000,
           0.1471,0.0400,0.0667,0.0930,0.2667,0.1538,0.2895,0.0811,0.1622,1.0000,
           0.0000,0.8333)



district <- factor(c(rep("DataII", length(data2)),
                     rep("DataIII", length(data3)),
                     rep("DataIV", length(data4))))

values <- c(data2, data3, data4)

# Kruskal-Wallis test
kruskal.test(values ~ district)

#install.packages("FSA")
library(FSA)
dunnTest(values ~ district, method="bonferroni")



library(ggplot2)
library(dplyr)
library(tidyr)

# Original data
interval_labels <- c("0.00-0.01","0.01-0.20","0.20-0.30","0.30-0.40","0.40-0.50",
                     "0.50-0.60","0.60-0.70","0.70-0.80","0.80-0.90","0.90-1.00")

#DataII  <- c(0.0010,0.0272,0.0241,0.0348,0.0502,0.0725,0.1048,0.1513,0.2185,0.3156)

DataII <- result$DataII
DataIII <- c(0.0041,0.0917,0.0606,0.0711,0.0833,0.0977,0.1146,0.1344,0.1576,0.1848)
DataIV  <- c(0.0397,0.5191,0.1510,0.1015,0.0682,0.0458,0.0308,0.0207,0.0139,0.0094)

# Combine into a tidy data frame
df <- data.frame(
  Interval = rep(interval_labels, times = 3),
  Probability = c(DataII, DataIII, DataIV),
  Dataset = rep(c("Data II", "Data III", "Data IV"), each = length(interval_labels))
)

windows(7,5)
# Plot with ggplot2
ggplot(df, aes(x = Interval, y = Probability, fill = Dataset)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(x = "Pass proportion interval", y = "Probability") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set1")


 # +coord_flip()




