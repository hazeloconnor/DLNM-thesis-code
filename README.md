# DLNM-thesis-code

Code repository for a thesis on Distributed Lag Non-Linear Models (DLNMs). Includes simulation studies in R using both frequentist and Bayesian approaches for parameter estimation.

## Structure

- `gaussian-model/`: Code for Gaussian DLNM simulations and estimation
- `nb-model/`: Code for Negative Binomial DLNM simulations, including a Polya-Gamma Gibbs sampler
- `PSI_FUN.dll`: Required compiled C library 

## External Code Notice

The `PSI_FUN.dll` file is compiled from source code by Daniel Dempsey ([Latent-Variable-DLM-Project](https://github.com/dempseydaniel/Latent-Variable-DLM-Project)) and is licensed under GNU GPL v3.0. Minor adjustments were made for compatibility. 
