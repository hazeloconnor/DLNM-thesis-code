## missing data
library(mice)
library(dlnm)
head (chicagoNMMAPS,1)
##first row of data
summary(chicagoNMMAPS)
##n/a tells us missing data, we see rhum has 1096 which is 21% do we leave it out? and pm10 has 251 only 5%
md.pattern(chicagoNMMAPS)
## tells us 3778 sets are complete
## 1085 are missing just rhum
## 240 are missing just pm10
## 11 are missing both
tempData<-mice(chicagoNMMAPS, m=5, maxit=50, meth='pmm', seed=500)
##m=5 is no of imputed data sets
## method is predictive mean matching
completeData<-complete(tempData,1)
##missing values replaced with imputed values