library(dlnm)
library(splines)

# Load complete data
source("impute_missing_data.R")

set.seed(123)


# Define parameters
B_0 <- 60            
variance <- 20     
num_sim <- 300          

# Define coefficients
theta1 <- c(0.2, -0.30, 0.10, -0.05)  
theta2 <- c(0.8, 0.2, -0.1,-0.1, 0.15)


# Cross-basis functions
cbpm10_1 <- crossbasis(completeData$pm10, lag = 25,
                       argvar = list(fun = "lin"),
                       arglag = list(fun = "poly", degree = 3))

cbpm10_2 <- crossbasis(completeData$pm10, lag = 25,
                       argvar = list(fun = "lin"),
                       arglag = list(fun = "ns", df = 4, knots=c(5,15,20)))


# Generate Linear Predictors
linear_predictor1 <- matrix(rep(B_0 + cbpm10_1 %*% theta1, num_sim), 
                            nrow = nrow(completeData), ncol = num_sim)
linear_predictor2 <- matrix(rep(B_0 + cbpm10_2 %*% theta2, num_sim), 
                            nrow = nrow(completeData), ncol = num_sim)


# Generate Error Terms
e_all1 <- matrix(rnorm(n = nrow(completeData) * num_sim, mean = 0, sd = sqrt(variance)),
                 nrow = nrow(completeData), ncol = num_sim)


# Simulated Response Variables
yi_1_all <- linear_predictor1 + e_all1
yi_2_all <- linear_predictor2 + e_all1




# Output:
# - yi_1_all: responses using polynomial lag structure
# - yi_2_all: responses using natural spline lag structure
# - cbpm10_1 and cbpm10_2: cross-basis matrices used in estimation
