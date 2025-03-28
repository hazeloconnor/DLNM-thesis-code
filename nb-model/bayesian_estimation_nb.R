library(BayesLogit)  
library(splines)
library(ggplot2)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(coda)


# This .dll is compiled from PSI_FUN.c (by Daniel Dempsey),
# Licensed under GNU GPL v3.0. This file was locally compiled for Windows compatibility.
dyn.load("PSI_FUN.dll") 
  

# Load data
source("nb-model/negative_binomial_simulation.R") 
Y <- simulated_Y 
data_sim <- data.frame(y = Y, Z [ ,-1])



##Function for priors
prior_function <- function( Z ) {
  p <- ncol(Z)
  prior <- list()
  
  # Prior for theta (regression coefficients)
  prior$theta$m <- rep(0, p)   
  prior$theta$v <- diag(100,p)     
  
  
  # Prior for r (overdispersion parameter)
  prior$r$a0 <- 5
  prior$r$b0 <- 1/10 
  
  #return result
  prior
  
}

## function for initial parameter for mcmc
initial_parameters <- function( Z ) {
  p <- ncol(Z)
  init <- list()
  
  # Initialize theta coefficients to zero
  init$theta <- rep(0, p)
  
  # Initialize overdispersion parameter r to 1
  init$r <- 1
  
  #return result
  init
}


# Gibbs Sampler Function
GibbsSampler_NB <- function(formula, data, prior, init, nsamples = 5000, burnin = 5000, thin = 1) {
  
  cat( "Initializing MCMC algorithm...\n" )
  total_iters <- nsamples + burnin
  Z_full <- model.matrix( formula , data)
  prior_r_a0 <- prior$r$a0
  prior_r_b0 <- prior$r$b0
  
  # Prior setup
  prior_cov_inv <- chol2inv(chol(prior$theta$v))  
  prior_precision_mean <- prior_cov_inv %*% prior$theta$m  
  
  
  # Storage
  p <- ncol( Z_full )
  theta_samples <- matrix(0, ncol = p, nrow = total_iters) 
  theta_samples[1, ] <- init$theta  
  r_samples <- numeric(total_iters)  
  r_samples[1] <- r <- init$r  
  
 
  colnames(theta_samples) <- colnames(Z_full)
  
  
  Z <- Z_full 
  Zb <- Z %*% init$theta 
  
  
  y <- data[[all.vars(formula)[1]]]
  n <- length( y )
  y_max <- max( y ) + 1
  
  
  # Initialize latent Polya-Gamma variables
  omega <- rpg( n , y + r, Zb )
  Zt0 <- t(Z * omega)
  lambda <- (y - r) / (2 * omega)
  V <- chol2inv(chol(prior_cov_inv + Zt0 %*% Z))
  B <- V %*% (prior_precision_mean + Zt0 %*% lambda)
  
  
  
  # Initialise matrix of pmfs for psi inference
  R <- matrix( 0, nrow = y_max, ncol = y_max )
  R[1, 1] <- R[2, 2] <- 1
  
  
  cat( 'Initialization complete. Running algorithm...\n' )
  for ( i in 2:total_iters ) {
    
    # Sample r
    sum_psi <- .C("PSI_FUN", 
                  as.double(r),             
                  as.integer(y),            
                  as.integer(n),       
                  as.integer(0:(y_max - 1)),
                  as.integer(y_max),        
                  as.double(R),             
                  as.integer(0)             
    )[[7]]  # Extract the output psi
    
    r_samples[i] <- r <- rgamma(1, 
                                prior_r_a0 + sum_psi, 
                                prior_r_b0 + sum(log(1 + exp(Zb))))
    
    
    # Sample theta
    omega <- rpg( n , y + r, Zb )
    
    Zt0 <- t(Z * omega)
    lambda <- (y - r) / (2 * omega)
    V <- chol2inv(chol(prior_cov_inv + Zt0 %*% Z))
    B <- V %*% (prior_precision_mean + Zt0 %*% lambda)
    theta_samples[i, ] <- B + t(chol(V)) %*% rnorm(ncol(Z))
    Zb <- Z %*% theta_samples[i, ]
    

    if (i %% 1000 == 0){
      cat( paste0( 'Completed iteration ', i, '...\n' ) )
    }
    
  } 
  
  cat('Algorithm complete. Returning result.\n')
  keep <- seq(burnin + 1, total_iters, thin)  
  list(
    theta = theta_samples[keep, ],  
    r = r_samples[keep],  
    init = init,  
    prior = prior,  
    Z = Z_full, 
    y = y  
  )
}

# Run the Sampler
set.seed(42)  
prior <- prior_function(Z)
init <- initial_parameters(Z)
nsamples <- 10000 
burnin  <- 5000  
thin <- 1  

results <- GibbsSampler_NB(y ~ ., data = data_sim, nsamples = nsamples, prior = prior, init = init, burnin = burnin, thin = thin)


# -------------------------------------------------
# Posterior Mean Estimates
# -------------------------------------------------

theta_posterior_mean <- colMeans(results$theta)
print(theta_posterior_mean)


r_posterior_mean <- mean(results$r)
print(r_posterior_mean)



# -------------------------------------------------
# MCMC Diagnostics: Trace Plots for r and intercept 
# -------------------------------------------------

# traceplot r
p1 <- ggplot(data = data.frame(iter = 1:length(results$r), r = results$r), aes(x = iter, y = r)) +
  geom_line(color = "blue") +
  ggtitle("Trace Plot of r") +
  xlab("Iteration") +
  ylab("r") +
  theme_minimal()

# traceplot intercept 
p2 <- ggplot(data = data.frame(iter = 1:nrow(results$theta), theta1 = results$theta[,1]), 
             aes(x = iter, y = theta1)) +
  geom_line(color = "blue") +
  ggtitle(expression("Trace Plot of Intercept " ~ theta[1])) +
  xlab("Iteration") +
  ylab(expression(theta[1])) +
  theme_minimal()

# Arrange plots in a 2-row, 1-column layout
grid.arrange(p1, p2, ncol = 1)



# Correlation between r and Intercept 
correlation_r_theta0 <- cor(results$r, results$theta[,1])
print(correlation_r_theta0)



# -------------------------------------------------
# Scatterplot: r vs. Intercept 
# -------------------------------------------------
plot_data <- data.frame(r = results$r, theta1 = results$theta[,1])

ggplot(plot_data, aes(x = r, y = theta1)) +
  geom_point(alpha = 0.5) +  # Transparency to show density
  geom_smooth(method = "lm", col = "red", se = FALSE) +  
  labs(title = "Scatterplot of r vs Intercept",
       x = "Overdispersion parameter (r)",
       y = expression(Intercept ~ (theta[1]))) +
  theme_minimal()



# -------------------------------------------------
# Distributed Lag Model: Estimating Lag Function
# -------------------------------------------------

lags <- 0:25
lag_basis_matrix <- ns(lags,knots = c(7,14), Boundary.knots = c(0, 25), intercept = TRUE)
true_theta_lag <-  theta[2:5]
true_lag_function <- lag_basis_matrix %*% true_theta_lag

num_samples <- 100  # Number of Posterior Samples to Plot
theta_samples <- results$theta[sample(1:nrow(results$theta), num_samples), 2:9]  

theta_samples_lag <- theta_samples[,1:4] 
lag_function_samples <- lag_basis_matrix %*% t(theta_samples_lag)

theta_mean <- theta_posterior_mean[2:5] 
lag_function_mean <- lag_basis_matrix %*% theta_mean

ylab_text <- expression(hat(gamma)[l])
plot(lags, lag_function_mean, type = "n", 
     ylim = range(lag_function_samples, true_lag_function), 
     xlab = "Lag (days)", ylab = ylab_text,
     main = "Estimated vs True Lag Function (Natural Spline Basis)")

for (i in 1:num_samples) {
  lines(lags, lag_function_samples[, i], col = rgb(0.7, 0.7, 0.7, 0.5))  # Light gray
}

lines(lags, lag_function_mean, col = "blue", lwd = 3)
lines(lags, true_lag_function, col = "red", lwd = 2, lty = 2)

legend("topright", legend = c("Posterior Samples", "Posterior Mean", "True Lag Function"), 
       col = c("grey", "blue", "red"), lwd = c(1, 3, 2), lty = c(1, 1, 2))
grid(lty = "dotted", col = "gray")


# -------------------------------------------------------------------
# MCMC Diagnostics (Trace + Density) remaining parameters
# -------------------------------------------------------------------
theta_samples <- as.mcmc(results$theta[, 2:9])
theta_df <- as.data.frame(as.matrix(theta_samples))


plot_trace_math <- function(theta_col, index) {
  df <- data.frame(Iteration = 1:length(theta_col), Value = theta_col)
  ggplot(df, aes(x = Iteration, y = Value)) +
    geom_line(color = "blue") +
    geom_hline(yintercept = mean(theta_col), color = "red", linetype = "dashed") +
    ggtitle(bquote("Traceplot of " ~ theta[.(index)])) +
    xlab("Iteration") + ylab("Value") +
    theme_minimal()
}

plot_density_math <- function(theta_col, index) {
  df <- data.frame(Value = theta_col)
  ggplot(df, aes(x = Value)) +
    geom_density(fill = "skyblue", alpha = 0.5, color = "blue") +
    geom_vline(xintercept = mean(theta_col), color = "red", linetype = "dashed") +
    ggtitle(bquote("Posterior Density of " ~ theta[.(index)])) +
    xlab("Value") + ylab("Density") +
    theme_minimal()
}


trace_plots <- lapply(1:8, function(i) plot_trace_math(theta_df[[i]], i + 1))
density_plots <- lapply(1:8, function(i) plot_density_math(theta_df[[i]], i + 1))


grid.arrange(grobs = trace_plots[1:4], ncol = 2)
grid.arrange(grobs = trace_plots[5:8], ncol = 2)

grid.arrange(grobs = density_plots[1:4], ncol = 2)
grid.arrange(grobs = density_plots[5:8], ncol = 2)

