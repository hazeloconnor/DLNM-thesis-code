library(dlnm)
library(splines)
library(MASS)
library(ggplot2)
library(gridExtra)


# Load data
source("gaussian-model/gaussian_simulated_data.R")
set.seed(123)


# Select outcome from first simulation
y <- yi_1_all[, 1]
X <- cbind(1, cbpm10_1)

complete_rows <- complete.cases(X, y)
X_clean <- X[complete_rows, ]
y_clean <- y[complete_rows]



#####################################
# 1. BAYESIAN ESTIMATION (GIBBS SAMPLING)
#####################################


# Prior hyperparameters
sigma2_theta <- 10     
alpha_prior <- 0.01    
beta_prior <- 0.01     
num_iter <- 10000      


# Initialize storage
n <- length(y_clean)
p <- ncol(X_clean)
theta <- matrix(0, nrow = num_iter, ncol = p) 
sigma2 <- numeric(num_iter)
theta[1, ] <- rep(0, p)  
sigma2[1] <- var(y_clean)


# Gibbs Sampler
for (i in 2:num_iter) {
  # Sample theta (coefficients)
  Q <- (1 / sigma2[i - 1]) * t(X_clean) %*% X_clean + 
    (1 / sigma2_theta) * diag(p)                         
  R <- chol(Q)                                           
  b <- (1 / sigma2[i - 1]) * t(X_clean) %*% y_clean      
  mu <- backsolve(R, backsolve(R, b, transpose = TRUE))  
  z <- rnorm(p)
  theta[i, ] <- mu + backsolve(R, z) 
  
  
  
  # Sample sigma^2 (variance)
  residuals <- y_clean - X_clean %*% theta[i, ]
  shape <- alpha_prior + (n / 2)
  rate <- beta_prior + (sum(residuals^2) / 2)
  sigma2[i] <- 1 / rgamma(1, shape, rate)
}


# Extract posterior samples
theta_samples <- theta[(num_iter / 2):num_iter, ]
sigma2_samples <- sigma2[(num_iter / 2):num_iter]

theta_mean <- colMeans(theta_samples)
sigma2_mean <- mean(sigma2_samples)


true_theta <- c(B_0, theta1)
squared_error_loss <- (true_theta - theta_mean)^2


cat("Posterior Means of Coefficients:\n")
print(theta_mean)
cat("\nPosterior Mean of Variance:\n")
print(sigma2_mean)
cat("\nSquared Error Loss (Posterior Mean vs. True Coefficients):\n")
print(squared_error_loss)


credible_intervals_theta <- apply(theta_samples, 2, quantile, probs = c(0.025, 0.975))
credible_interval_sigma2 <- quantile(sigma2_samples, probs = c(0.025, 0.975))

cat("\n95% Credible Intervals for Theta:\n")
for (j in 1:ncol(theta_samples)) {
  cat(sprintf("Theta_%d: Mean = %.5f, 95%% CI = [%.5f, %.5f]\n",
              j, theta_mean[j], credible_intervals_theta[1, j], credible_intervals_theta[2, j]))
}

cat(sprintf("\nSigma^2: Mean = %.5f, 95%% CI = [%.5f, %.5f]\n",
            sigma2_mean, credible_interval_sigma2[1], credible_interval_sigma2[2]))



#####################################
# 2. CONVERGENCE DIAGNOSTICS
#####################################

# Trace plots
theta_plots <- list()

for (j in 1:ncol(theta_samples)) {
  theta_mean_j <- mean(theta_samples[, j])  
  
  p <- ggplot(data.frame(iter = 1:nrow(theta_samples), value = theta_samples[, j]),
              aes(x = iter, y = value)) +
    geom_line(color = "blue") +
    geom_hline(yintercept = theta_mean_j, color = "red", linetype = "dashed") +
    ggtitle(bquote("Trace Plot for " ~ theta[.(j)])) +  
    xlab("Iteration") + 
    ylab("Value") +  
    theme_minimal()
  
  theta_plots[[j]] <- p  
}

# Trace plot for sigma^2
sigma2_plot <- ggplot(data.frame(iter = 1:length(sigma2_samples), value = sigma2_samples),
                      aes(x = iter, y = value)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = sigma2_mean, color = "red", linetype = "dashed") +  
  ggtitle(expression("Trace Plot for " ~ sigma^2)) + 
  xlab("Iteration") + 
  ylab("Value") + 
  theme_minimal()

theta_plots[[length(theta_plots) + 1]] <- sigma2_plot

grid.arrange(grobs = theta_plots, ncol = 2)



# Density plots
plots <- list()
plot_density <- function(theta_col, index) {
  df <- data.frame(Value = theta_col)
  ggplot(df, aes(x = Value)) +
    geom_density(fill = "skyblue", alpha = 0.5, color = "blue") +
    geom_vline(xintercept = mean(theta_col), color = "red", linetype = "dashed", lwd = 1) +
    ggtitle(bquote("Posterior Density of " ~ theta[.(index)])) +  
    xlab("Value") +
    ylab("Density") +
    theme_minimal()
}

for (j in 1:ncol(theta_samples)) {
  p <- plot_density(theta_samples[, j], j)  # Pass j to function
  plots[[j]] <- p
}


# Density plot for sigma^2
sigma2_hist <- ggplot(data = data.frame(value = sigma2_samples), aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.5, color = "blue") +
  ggtitle(expression("Posterior Density of " ~ sigma^2)) +
  xlab("Value") + 
  ylab("Density") +
  geom_vline(xintercept = mean(sigma2_samples), color = "red", linetype = "dashed", lwd = 1) +
  theme_minimal()

plots[[length(plots) + 1]] <- sigma2_hist
grid.arrange(grobs = plots, ncol = 2)



#####################################
# 3. PLOTTING LAG FUNCTION
#####################################

lags <- 0:25
scaled_lags <- lags / 25  
lag_basis_poly <- poly(scaled_lags, degree = 3, raw = TRUE) 
lag_basis_poly <- cbind(Intercept = 1, lag_basis_poly)


true_lag_function <- lag_basis_poly %*% theta1 
num_samples <- 200  
theta_samples <- theta[sample((num_iter / 2):num_iter, num_samples), 2:5] 
lag_function_samples <- lag_basis_poly %*% t(theta_samples) 
lag_function_mean <- lag_basis_poly %*% theta_mean[2:5]


ylab_text <- expression(hat(gamma)[l])

plot(lags, lag_function_mean, type = "n", 
     ylim = range(lag_function_samples, true_lag_function), 
     xlab = "Lag (days)", 
     ylab = ylab_text, 
     main = "Estimated vs True Lag Function (Polynomial Basis)")

for (i in 1:num_samples) {
  lines(lags, lag_function_samples[, i], col = rgb(0.7, 0.7, 0.7, 0.5))
}

lines(lags, lag_function_mean, col = "blue", lwd = 3)
lines(lags, true_lag_function, col = "red", lwd = 2, lty = 2)
abline(h = 0, lty = 2, col = "black")
grid(lty = "dotted", col = "gray")

legend("topright", 
       legend = c("True Lag Function", 
                  "Posterior Mean", 
                  "Posterior Samples"),
       col = c("red", "blue", rgb(0.7, 0.7, 0.7)), 
       lty = c(2, 1, 1), 
       lwd = c(2, 3, 1), 
       bty = "o",  # box legend
       cex = 0.9)

