
library(mice)
library(lavaan)
library(dplyr)
library(readr)



# (2) load data -----------------------------------------------------------
complete_data <- read_csv("YEPPersonality_all data_122423_nocarelessblv2.csv")

select_data <- complete_data %>% 
  select(ddad_csrtot_cur_t5, consc_t5, open_t5)




# evaluate the mean and missing rate of the data --------------------------
round(mean(select_data$ddad_csrtot_cur_t5, na.rm = T),2)
round(mean(select_data$consc_t5, na.rm = T),2)
round(mean(select_data$open_t5, na.rm = T),2)

round(sum(is.na(select_data$ddad_csrtot_cur_t5))/nrow(select_data),2)
round(sum(is.na(select_data$consc_t5))/nrow(select_data),2)
round(sum(is.na(select_data$open_t5))/nrow(select_data),2)

# (4) apply fiml and mice -------------------------------------------------

# fiml --------------------------------------------------------------------


model_fiml <- try(sem(model = "ddad_csrtot_cur_t5~ consc_t5 + open_t5", 
                      data = select_data, 
                      missing = "fiml", 
                      fixed.x = FALSE))

(result_fiml <- summary(model_fiml))




# listwise ----------------------------------------------------------------

miss_indicator <- numeric(length = nrow(select_data))

for (i in 1:nrow(select_data)) {
  miss_indicator[i] <- sum(is.na(select_data[i, ]))
}

# Find rows with no missing values
complete_index <- which(miss_indicator == 0)

complete_data <- try(select_data[complete_index, ])

model <- lm(ddad_csrtot_cur_t5~ consc_t5 + open_t5 , data = complete_data)

summary_model <- summary(model)

summary_model$coefficients

# cart --------------------------------------------------------------------


complete_data_cart <- try(mice(data = select_data,
                               m = 5, # generate 5 imputed data
                               method = "cart",
                               maxit = 5,
                               seed = 20251102,
                               printFlag = F))


(result_cart <- summary(pool(with(complete_data_cart, lm(ddad_csrtot_cur_t5~ consc_t5 + open_t5)))))




# rf ----------------------------------------------------------------------
complete_data_rf <- try(mice(data = select_data,
                             m = 5, # generate 5 imputed data
                             method = "rf",
                             maxit = 5,
                             seed = 20251102,
                             printFlag = F))


(result_rf <- summary(pool(with(complete_data_rf, lm(ddad_csrtot_cur_t5~ consc_t5 + open_t5)))))



# pmm ---------------------------------------------------------------------

complete_data_pmm <- try(mice(data = select_data,
                              m = 5, # generate 5 imputed data
                              method = "pmm",
                              maxit = 5,
                              seed = 20251102, 
                              printFlag = F))


(result_pmm <- summary(pool(with(complete_data_pmm, lm(ddad_csrtot_cur_t5~ consc_t5 + open_t5)))))



# norm --------------------------------------------------------------------

complete_data_norm <- try(mice(data = select_data,
                               m = 5, # generate 5 imputed data
                               method = "norm",
                               maxit = 5,
                               seed = 20251102,
                               printFlag = F))


(result_norm <- summary(pool(with(complete_data_norm, lm(ddad_csrtot_cur_t5~ consc_t5 + open_t5)))))




# (5) calculate the estimate for simulation  ------------------------------

coef_1 <- round(mean(result_fiml$pe[1,5],summary_model$coefficients[2,1],result_cart[2,2],result_rf[2,2],result_pmm[2,2],result_norm[2,2]),3)



coef_2 <- round(mean(result_fiml$pe[2,5],summary_model$coefficients[3,1],result_cart[3,2],result_rf[3,2],result_pmm[3,2],result_norm[3,2]),3)



# (6) simulate the data ---------------------------------------------------

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
