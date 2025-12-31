library(lavaan)
library(psych)
library(dplyr)
library(mice)

# T1 five charc
dataT115 <- read.csv('T1 T15 Big5_YEPResSample2.csv',header = T, na.strings = c("*", "NA"))
colnames(dataT115) <- c("id","e11","e12","e13","e14","e15","e16","e17","e18", "a11","a12","a13","a14","a15", 
                        "a16","a17","a18","c11","c12","c13","c14","c15","c16","c17", "n11", "n12", "n13", "n14",
                        "n15","n16", "n17","n18", "o11", "o12", "o13", "o14", "o15", "o16", "o17", "o18", 
                        "e81","e82","e83","e84","e85","e86","e87","e88", "a81","a82","a83","a84","a85", 
                        "a86","a87","a88","c81","c82","c83","c84","c85","c86","c87", "n81", "n82", "n83", "n84",
                        "n85","n86", "n87","n88", "o81", "o82", "o83", "o84", "o85", "o86", "o87", "o88")

dfec <- dataT115[, c("e11","e12","e13","e14","e15","e16","e17","e18", 
                     "c81","c82","c83","c84","c85", "c86","c87")]

#missing rate
proportion <- round(100*apply(dfec, 2, function(x) mean(is.na(x))),3)
proportion

describe(dfec)

#EC
SEM_ec <- 'F1  =~ e11 + e12 + e13 + e14 + e15 + e16 + e17 + e18
                F2  =~ c81 + c82 + c83 + c84 + c85 + c86 + c87
                F2 =~ F1'

sem_ec <- sem(SEM_ec, data = dataT115, std.lv = TRUE,missing = "fiml")
summary(sem_ec,fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
s_sem_ec <- standardizedSolution(sem_ec)
s_sem_ec

sem_ec_list <- sem(SEM_ec, data = dataT115, std.lv = TRUE,missing = "listwise")
summary(sem_ec_list,fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
s_sem_ec_list <- standardizedSolution(sem_ec_list)
s_sem_ec_list
