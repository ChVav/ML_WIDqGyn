### testing August 15 - output App ----
# only simplest case implemented (pmr calc with on-plate cal curve)

library(dplyr)
library(tidyverse)
library(openxlsx)

# read in output app
unzip("../../TestDataqWidGyn/HALL/3-outputTEST/Final_results_2022-08-13.zip", exdir="../../TestDataqWidGyn/HALL/3-outputTEST/Final_results_2022-08-13/")
list.files("../../TestDataqWidGyn/HALL/3-outputTEST/Final_results_2022-08-13/")
test <- read.xlsx("../../TestDataqWidGyn/HALL/3-outputTEST/Final_results_2022-08-13/batch_results.xlsx", sheet=1)
test <- test %>% filter(!Sample %in% c("STD1", "STD2", "STD3", "STD4", "gBlock", "H2O", "STD_1", "STD_2", "STD_3", "STD_4", "gBLOCK"))

# read in expected data
data_expected <- read.xlsx("../../TestDataqWidGyn/HALL/3-outputCHECK2/checkOnplateCalib.xlsx", sheet=1)
qEC_threshold1 <- 0.03
qEC_threshold2 <- 0.63
data_expected <- data_expected %>% #42 samples, no GYPC2 in samples amplified
  mutate(qEC = (GYPC1 + GYPC2 + ZSCAN12),
         qCIN = (RALYL + DPP6 + GSX1)) %>%
  filter(!Sample %in% c("STD1", "STD2", "STD3", "STD4", "gBlock", "H2O", "STD_1", "STD_2", "STD_3", "STD_4", "gBLOCK")) #42 samples GYPC2 was NA?
data_expected <- data_expected %>%
  mutate(result_qEC = case_when(data_expected$qEC < qEC_threshold1 ~ "low risk EC/CIN",
                                data_expected$qEC >= qEC_threshold1 & data_expected$qEC <= qEC_threshold2 ~ "high risk EC/CIN",
                                data_expected$qEC > qEC_threshold2 ~ "very high risk EC/CIN"))

# compare results
identical(test$Sample,data_expected$Sample) #TRUE
identical(test$WIDqEC_test, data_expected$result_qEC) #TRUE :)
identical(round(test$WIDqEC,10), round(data_expected$qEC,10)) #TRUE, assume differences in rounding
identical(round(test$WIDqCIN,10), round(data_expected$qCIN,10)) #TRUE, assume differences in rounding

### testing August 11 - main function: run calculate_pmr_v3.R separately and compare with output pipeline Chiara ----
# corrected multipoint regression
# corrected means
# use only data as input
# add sum PMRs and result qEC test based on thresholds
# plate-calibration curve is downloaded as well

library(dplyr)
library(tidyverse)
library(openxlsx)
library(ggplot2)
library(ggpubr)

source("./calculate_pmr_v3.R")

data <- read.table(file= "../../TestDataqWidGyn/HALL/2-outputCHECK2/test_Results_20220808 141627.csv",
                    sep = ",",
                    header = TRUE) # data exported with ALL thresholds set manually
#curve <- read.table(file= "../../TestDataqWidGyn/HALL/2-outputCHECK2/test_Standard Curve Result_20220808 141627.csv",
#                   sep = ",",
#                   header = TRUE) # data exported with ALL thresholds set manually

####### On plate-regression curve
### for test data calculate PMR, the WID-qEC and WID-qCIN, and outcome for WID-qEC ----

test <- calculate_pmr(data)
test <- test[[2]] %>% filter(!Sample %in% c("STD1", "STD2", "STD3", "STD4", "gBlock", "H2O", "STD_1", "STD_2", "STD_3", "STD_4", "gBLOCK")) #42 samples no GYPC2, no GYPC2 in samples amplified

### Comparison expected vs current output for on-plate calibration---- 
# note, visual inspection calibration curve = OK
data_expected <- read.xlsx("../../TestDataqWidGyn/HALL/3-outputCHECK2/checkOnplateCalib.xlsx", sheet=1)
data_expected <- data_expected %>% #42 samples, no GYPC2 in samples amplified
  mutate(qEC = (GYPC1 + GYPC2 + ZSCAN12),
         qCIN = (RALYL + DPP6 + GSX1)) %>%
  filter(!Sample %in% c("STD1", "STD2", "STD3", "STD4", "gBlock", "H2O", "STD_1", "STD_2", "STD_3", "STD_4", "gBLOCK")) #42 samples GYPC2 was NA?
data_expected <- data_expected %>%
  mutate(result_qEC = case_when(data_expected$qEC < qEC_threshold1 ~ "low risk EC/CIN",
                                data_expected$qEC >= qEC_threshold1 & data_expected$qEC <= qEC_threshold2 ~ "high risk EC/CIN",
                                data_expected$qEC > qEC_threshold2 ~ "very high risk EC/CIN"))
identical(test$Sample,data_expected$Sample) #TRUE
identical(test$WIDqEC_test, data_expected$result_qEC) #TRUE :)
identical(round(test$WIDqEC,10), round(data_expected$qEC,10)) #TRUE, assume differences in rounding
identical(round(test$WIDqCIN,10), round(data_expected$qCIN,10)) #TRUE, assume differences in rounding

####### Fixed regression coeff
test2 <- calculate_pmr(data,calib_fixed=TRUE)
test2 <- test2[[2]] %>% filter(!Sample %in% c("STD1", "STD2", "STD3", "STD4", "gBlock", "H2O", "STD_1", "STD_2", "STD_3", "STD_4", "gBLOCK")) #42 samples, no GYPC2 in samples amplified

# note: using on-plate calibration (R2=0.994) or fixed regression coefficients does not give the same results!!!!!
identical(test$Sample,test2$Sample) #TRUE
identical(test$WIDqEC_test, test2$WIDqEC_test) #FALSE

#3522103497
#3522103435
#3522103426
#3522103468

### Comparison expected vs current output for fixed regression coefficients ---- 
data_expected2 <- read.xlsx("../../TestDataqWidGyn/HALL/3-outputCHECK2/checkFixedCoeff.xlsx", sheet=1)
data_expected2 <- data_expected2 %>% #42 samples, no GYPC2 in samples amplified
  mutate(qEC = (GYPC1 + GYPC2 + ZSCAN12),
         qCIN = (RALYL + DPP6 + GSX1)) %>%
  filter(!Sample %in% c("STD1", "STD2", "STD3", "STD4", "gBlock", "H2O", "STD_1", "STD_2", "STD_3", "STD_4", "gBLOCK")) #42 samples GYPC2 was NA?
data_expected2 <- data_expected2 %>%
  mutate(result_qEC = case_when(data_expected2$qEC < qEC_threshold1 ~ "low risk EC/CIN",
                                data_expected2$qEC >= qEC_threshold1 & data_expected2$qEC <= qEC_threshold2 ~ "high risk EC/CIN",
                                data_expected2$qEC > qEC_threshold2 ~ "very high risk EC/CIN"))
identical(test2$Sample,data_expected2$Sample) #TRUE
identical(test2$WIDqEC_test, data_expected2$result_qEC) #TRUE :)
identical(test2$WIDqEC, data_expected2$qEC) #FALSE, small differences cannot be due to data export...

identical(round(test2$WIDqEC,10), round(data_expected2$qEC,10)) #TRUE, assume differences in rounding
identical(round(test2$WIDqCIN,10), round(data_expected2$qCIN,10)) #TRUE, assume differences in rounding

