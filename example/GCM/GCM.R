library(lavaan)
library(psych)
library(dplyr)
library(mice)

#all data
dataall <- read.csv('YEPPersonality_all data_122423_nocarelessblv3.csv',header = T, na.strings = c(" ", "NA"))

data <- dataall[, c("id","Hollingshead_SES_t1", "age", "gender", "ethnicity_categorical", "gnf.1", "gnf.3", 
                    "gnf.5", "gnf.7", "neur_t1", "neur_t3", "neur_t5", "neur_t7", "extra_t1", "extra_t3", 
                    "extra_t5", "extra_t7", "consc_t1", "consc_t3", "consc_t5", "consc_t7", "agree_t1", 
                    "agree_t3", "agree_t5", "agree_t7", "open_t1", "open_t3", "open_t5", "open_t7",
                    "cs_sum_t1", "cs_sum_t3", "cs_sum_t5", "cs_sum_t7", "ddad_csrtot_cur_t1",
                    "ddad_csrtot_cur_t3", "ddad_csrtot_cur_t5", "ddad_csrtot_cur_t7", "Zcs_sum_t1",           
                    "Zcs_sum_t3", "Zcs_sum_t5", "Zcs_sum_t7", "Zddad_csrtot_cur_t1", "ddad_cp_z3", "ddad_cp_z5", "ddad_cp_z7")]


#agree
dataagree <- data[, c("agree_t1", "agree_t3", "agree_t5", "agree_t7")]
colnames(dataagree) <- c("V1","V2","V3","V4")

# descriptive statistics
describe(dataagree)

#missing rate
proportion <- round(100*apply(dataagree, 2, function(x) mean(is.na(x))),3)
proportion

gcmodel<-'i =~ 1*V1 + 1*V2 + 1*V3  + 1*V4
	s =~ 0*V1 + 1*V2 + 2*V3 + 3*V4'

res.agree <- growth(gcmodel, dataagree, missing = "fiml")
parameterEstimates(res.agree)[c(13:15,20:21),]

#listwise
res.agree2 <- growth(gcmodel, dataagree, missing = "listwise")
parameterEstimates(res.agree2)[c(13:15,20:21),]




