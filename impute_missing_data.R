## Imputation of missing values in the chicagoNMMAPS dataset 
library(mice)
library(dlnm)

head (chicagoNMMAPS,1)
summary(chicagoNMMAPS)

# Check missing data structure
md.pattern(chicagoNMMAPS)

# Notes:
# - 'rhum' has 1,096 missing values (~21%)
# - 'pm10' has 251 missing (~5%)
# Decision: Use multiple imputation to retain rhum and pm10.

# Perform multiple imputation (5 datasets)
tempData<-mice(chicagoNMMAPS, m=5, maxit=50, meth='pmm', seed=500)

# Extract the first complete dataset
completeData<-complete(tempData,1)
##missing values replaced with imputed values