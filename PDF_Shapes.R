
# Define the PDF function
pdf_theta <- function(x, theta) {
  if (theta == 1) stop("theta cannot be exactly 1")
  (log(theta) / (theta - 1)) * (theta^x)
}

# Grid of x values
x_vals <- seq(0, 1, length.out = 1000)

# Parameter values
theta_values <- c(0.005, 0.05, 1.05, 5, 50, 500)

# Empty plot
windows(5, 5)
par(mar = c(4.5, 4.5, 1.1, 1.1))
plot(x_vals, pdf_theta(x_vals, theta_values[1]), type = "n",
     xlab = "x",
     ylab = expression(f(x*"; "*theta)),
     ylim = c(0, 5))

# Add curves
colors <- c("red","blue","darkgreen","purple",
  "orange","brown")
#colors <- rainbow(length(theta_values))
for (i in seq_along(theta_values)) {
  theta <- theta_values[i]
  y_vals <- pdf_theta(x_vals, theta)
  lines(x_vals, y_vals, col = colors[i], lwd = 2)
}

# Legend with Greek theta
legend_labels <- sapply(theta_values, function(th)
  bquote(theta == .(th))
)

legend("top", cex = 0.7,
       legend = legend_labels,
       col = colors,
       lwd = 2,
       bty = "T")

#-------------------------------------##############
######## PDF Plots #################################
#------------------------------------- 

###################################################
# Extended PDF with parameter theta (>0, theta != 1)
pdf_theta <- function(x, theta) {
  if (theta == 1) stop("theta cannot be exactly 1")
  (log(theta) / (theta - 1)) * (theta^x)
}
  # Grid of x values
  x_vals <- seq(0, 1, length.out = 1000)
  
  # Define a set of parameter combinations to illustrate flexibility
  # alpha = 0.5, 1, 2 ; theta = 0.2, 0.8, 1.2, 5
  param_combos <- c(0.005, 0.05, 1.05, 5, 50, 500)
  
  windows(6,4)
  # Set up a multi-panel plot (3 rows, 4 columns)
  par(mfrow = c(2, 3), mar = c(2, 2, 2, 1), oma = c(2, 2, 0, 0))
  
  for (i in 1:length(param_combos)) {
    theta <- param_combos[i]
    
    # Compute density
    y_vals <- pdf_theta(x_vals, theta)
    
    # Plot
    plot(x_vals, y_vals, type = "l", lwd = 2, col = "steelblue",
         xlab = "", ylab = "", 
         main = bquote(theta == .(theta)))
    grid()
  }
  
  # Add overall labels
  mtext("x", side = 1, outer = TRUE, line = 0.5)
  mtext(expression(f(x*"; "*theta)), side = 2, outer = TRUE, line = 0.3)
  
#-------------------------------------##############
######## HRF Plots #################################
#------------------------------------- 
hrf <- function(x, theta) {
    (log(theta) * theta^x) / (theta - theta^x)
  }
  # Grid of x values
  x_vals <- seq(0, 1, length.out = 1000)
  
  # Define a set of parameter combinations to illustrate flexibility
  # alpha = 0.5, 1, 2 ; theta = 0.2, 0.8, 1.2, 5
  param_combos <- c(0.005, 0.05, 1.05, 5, 50, 500)
  
   windows(6,4)
  # Set up a multi-panel plot (3 rows, 4 columns)
  par(mfrow = c(2, 3), mar = c(2, 2, 2, 1), oma = c(2, 2, 0, 0))
  
  for (i in 1:length(param_combos)) {
    theta <- param_combos[i]
    
    # Compute density
    y_vals <- hrf(x_vals, theta)
    
    # Plot
    plot(x_vals, y_vals, type = "l", lwd = 2, col = "blue",
         xlab = "", ylab = "", 
         main = bquote(theta == .(theta)))
    grid()
  }
  
  # Add overall labels
  mtext("x", side = 1, outer = TRUE, line = 0.5)
  mtext(expression(h(x*";"*theta)), side = 2, outer = TRUE, line = 0.3)
  
  


