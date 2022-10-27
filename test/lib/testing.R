
#testing main workfunctions

library(dplyr)
library(tidyverse)
library(ggpubr)
library(openxlsx)
library(ggplot2)

#### cal curve ----
source("../../lib/calibGblock.R")

calib_fixed <- FALSE

# read in dummy data
data <- read.table(file = "../../data/dummy/1-input/test_Results_20220808 141627.csv",
                   sep = ",",
                   header = TRUE)
calibration <- calibGblock(data)
calibration[[1]] #gives plot
intercept <- calibration[[2]] #gives intercept
slope <- calibration[[3]] # gives slope

#### calcPMRGyn.R ----
source("../../lib/calcPMRGyn.R")

# collect target and sample overviews
targets <- unique(data$Target)
samples <- unique(data$Sample)

data <- data %>%
  mutate(rep = case_when(grepl("A|C|E|G|I|K|M|O", Well.Position) == TRUE ~ 1,
                         grepl("B|D|F|H|J|L|N|P", Well.Position) == TRUE ~ 2))

threshold_COL2A1=30
threshold_targets=38

calcPMR1 <- calcPMRGyn(data,targets,samples,intercept, slope,threshold_COL2A1, threshold_targets)

#### calcPMRGyn_extcurve.R ----
#do I get the same results when using same plate as ext_curve
source("../../lib/calcPMRGyn_extcurve.R")

samplesdata <- data
gblockdata <- read.table(file = "../../data/dummy/1-input/test_Standard Curve Result_20220808 141627.csv",
                         sep = ",",
                         header = TRUE)

calcPMR2 <- calcPMRGynExt(samplesdata,gblockdata,intercept, slope,threshold_COL2A1, threshold_targets)
identical(calcPMR1[[1]],calcPMR2[[1]]) #TRUE!!!!!!!!!!!!!! Halleluja :)
