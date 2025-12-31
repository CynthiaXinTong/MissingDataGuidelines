
library(mice)
library(lavaan)
library(dplyr)
library(readr)


generate_data <- function(N, rep) {
  set.seed(rep)
  x1 <- rnorm(N)
  x2 <- rnorm(N)
  y <- -0.31 * x1 + 0.141 * x2 + rnorm(N) # coef_1=-0.31   coef_2=0.141
  
  data <- as.data.frame(cbind(y, x1, x2)) # depend on how many auxiliary variables need, extract from auxiliary variable data frame
  return(data)
}

data <- generate_data(N = 100, rep = 114514)

summary(lm(y~x1+x2, data=data))
