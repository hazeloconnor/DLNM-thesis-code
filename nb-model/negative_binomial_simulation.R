library(dlnm)
library(MASS)

# Load data
source("impute_missing_data.R")
set.seed(123)


r <- 50  # Dispersion parameter

theta <- c(-1,   # Intercept
           2.5, -1.5, 1.5, -0.7,  
           1.5, -1, 2, -0.5)  


# Construct crossbasis
cbpm10_3 <- crossbasis(
  completeData$pm10, 
  lag = 25,
  argvar = list(fun = "poly", degree=2),
  arglag = list(fun = "ns", df = 4, knots = c(7, 14))
)

# Construct design matrix
Z <- cbind(1, as.matrix(cbpm10_3))

# Compute linear predictor and filter valid rows
log_mu <- Z %*% theta
mu <- exp(log_mu)
valid_rows <- !is.na(mu)
mu <- mu[valid_rows]
Z <- Z[valid_rows, ]



# Simulate Negative Binomial response
simulated_Y <- rnbinom(n = length(mu), size = r, mu = mu)


mean_simulated_Y <- mean(simulated_Y)
var_simulated_Y <- var(simulated_Y)

cat("Mean of simulated Y:", mean_simulated_Y, "\n")
cat("Variance of simulated Y:", var_simulated_Y, "\n")

