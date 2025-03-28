library(dlnm)
library(splines)
library(MASS)

# Load data
source("gaussian-model/gaussian_simulated_data.R")
set.seed(123)



#####################################
# 1. LEAST SQUARES ESTIMATION (LSE)
#####################################


# Construct Design Matrices
X1 <- cbind(1, cbpm10_1)  
X2 <- cbind(1, cbpm10_2)  


# Least Squares Estimation via Cholesky decomposition
estimate_lse_backsolve <- function(X, Y) {
  complete_cases <- complete.cases(X, Y)
  X_clean <- X[complete_cases, , drop = FALSE]
  Y_clean <- Y[complete_cases, ]
  
  # Compute Cholesky decomposition
  XtX <- t(X_clean) %*% X_clean
  chol_XtX <- chol(XtX)
  
  # Solve for theta using backsolve
  theta_hat <- backsolve(chol_XtX, backsolve(chol_XtX, t(X_clean) %*% Y_clean, transpose = TRUE))
  
  return(theta_hat)
}


theta_hats_1 <- estimate_lse_backsolve(X1, yi_1_all)
theta_hats_2 <- estimate_lse_backsolve(X2, yi_2_all)


# Compute Mean and Standard Deviation of Estimates
theta_mean_1 <- rowMeans(theta_hats_1)
theta_sd_1 <- apply(theta_hats_1, 1, sd)

theta_mean_2 <- rowMeans(theta_hats_2)
theta_sd_2 <- apply(theta_hats_2, 1, sd)

# Print Results
cat("Model 1: Coefficients (Polynomial Lag Basis)\n")
cat("Mean coefficients:\n"); print(theta_mean_1)
cat("Standard deviation of coefficients:\n"); print(theta_sd_1)

cat("\nModel 2: Coefficients (Natural Splines Lag Basis)\n")
cat("Mean coefficients:\n"); print(theta_mean_2)
cat("Standard deviation of coefficients:\n"); print(theta_sd_2)



theta1_true <- c(60, 0.2, -0.30, 0.10, -0.05)
theta2_true <- c(60, 0.8, 0.2, -0.1, -0.1, 0.15)

# Squared Error Loss Function Q 
Q1 <- (theta1_true - theta_mean_1)^2
Q2 <- (theta2_true - theta_mean_2)^2


cat("\nModel 1: Squared Error Loss Q\n"); print(Q1)
cat("\nModel 2: Squared Error Loss Q\n"); print(Q2)


#####################################
# 2. PREDICTION OF EFFECTS
#####################################


lags <- 0:25  

## 2.1 Prediction Using Polynomial Basis (Model 1)
scaled_lags <- lags / 25  
lag_basis_poly <- poly(scaled_lags, degree = 3, raw = TRUE)  
lag_basis_poly <- cbind(1, lag_basis_poly)  # Include intercept
lag_effects_poly <- lag_basis_poly %*% theta_mean_1[-1]  

# Compute Variance and Confidence Intervals
var_theta_poly <- cov(t(theta_hats_1))  
var_lag_effects_poly <- lag_basis_poly %*% var_theta_poly[-1, -1] %*% t(lag_basis_poly)
se_lag_effects_poly <- sqrt(diag(var_lag_effects_poly))
lower_ci_poly <- lag_effects_poly - 1.96 * se_lag_effects_poly
upper_ci_poly <- lag_effects_poly + 1.96 * se_lag_effects_poly

# Plot Prediction for Polynomial Basis
par(mar = c(5, 6, 4, 2))
plot(lags, lag_effects_poly, type = "l", lwd = 3, col = "blue",
     xlab = "Lag (Days)", ylab = expression(hat(gamma)),
     main = "Polynomial Lag-Response Function with 95% CI")
polygon(c(lags, rev(lags)), c(lower_ci_poly, rev(upper_ci_poly)), col = adjustcolor("blue", alpha.f = 0.2), border = NA)
lines(lags, lower_ci_poly, col = "red", lty = 2, lwd = 2)
lines(lags, upper_ci_poly, col = "red", lty = 2, lwd = 2)
abline(h = 0, lty = 2, col = "black")
grid(lty = "dotted", col = "gray")




## 2.2 Prediction Using Natural Splines (Model 2)
lag_basis_ns <- ns(lags, df = 4, Boundary.knots = c(0, 25), knots = c(5, 15, 20), intercept = TRUE)
lag_effects_ns <- lag_basis_ns %*% theta_mean_2[-1]

# Compute Variance and Confidence Intervals
var_theta_ns <- cov(t(theta_hats_2))  
var_lag_effects_ns <- lag_basis_ns %*% var_theta_ns[-1, -1] %*% t(lag_basis_ns)
se_lag_effects_ns <- sqrt(diag(var_lag_effects_ns))
lower_ci_ns <- lag_effects_ns - 1.96 * se_lag_effects_ns
upper_ci_ns <- lag_effects_ns + 1.96 * se_lag_effects_ns

# Plot Prediction for Natural Splines
par(mar = c(5, 6, 4, 2))
plot(lags, lag_effects_ns, type = "l", lwd = 3, col = "blue",
     xlab = "Lag (Days)", ylab = expression(hat(gamma)),
     main = "Natural Splines Lag-Response Function with 95% CI")
polygon(c(lags, rev(lags)), c(lower_ci_ns, rev(upper_ci_ns)), col = adjustcolor("blue", alpha.f = 0.2), border = NA)
lines(lags, lower_ci_ns, col = "red", lty = 2, lwd = 2)
lines(lags, upper_ci_ns, col = "red", lty = 2, lwd = 2)
abline(h = 0, lty = 2, col = "black")
grid(lty = "dotted", col = "gray")

  
  

