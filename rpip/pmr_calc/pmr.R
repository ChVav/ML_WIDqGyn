# PMR calculation
# Author: Chiara Herzog
# contact: chiara.herzog@uibk.ac.at
# Date 17 Nov 2021
# Update log:
#   - CV 08 Aug 2022: Fix so that output folder is created when needed.
#   - CH 26 Jul 2022: Amend pipeline so one PMR can be carried forward
#   - CH 25 Jul 2022: Paths can either be relative (default - in raw data folder) or FIXED (input from other locations possible)
#   - CV 06 Jul 2022: Added sheet to output with estimated input DNA
#   - CV 30 Jun 2022: Version for incomplete WID-qGYN, option for using fixed regression coefficients and COL2A1 Ct thresholds
# Description:
#   This version works for WID-qGYN, with minimally COL2A1, EpC and ImC targets. Targets to exclude
#   are set in the previous part of the pipeline (sheet generation).
#   This is a wrapper function. It calls the most recent PMR script and calculates PMRs as well as
#   performs an autotest for the most up to date script to see whether it replicates previous
#   version's values on a barcelona plate.
# Variables:
#   - calib_fixed: Instead of using the on-plate calibration curve, fixed values for intercept and slope
#     can be used. To use a previous (single) calibration curve, replace exported calibration file from
#     a different plate.
#   - threshold_COL2A1: Optionally, a threshold can be set for COL2A1, default = 35 (everything higher fails).
#     For all other targets, Ct is set to NA when Amp status is inconclusive, but one Ct/PMR per rep is
#     carried forward.
#   - type: gblock or sss1, but sss1 has not been implemented yet.
#   - path: relative (default: location in raw data folder) or fixed (input from other locations)


pmr <- function(folder, output, experimentname,
                write.results = TRUE,
                type = "gblock",
                calib_fixed = FALSE,
                fix_intercept = 36.9,
                fix_slope = -3.4,
                threshold_COL2A1 = 35,
                path = "relative"){
  
  # load libraries
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(openxlsx))
  suppressPackageStartupMessages(library(fs))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggpubr))
  
  # specify paths input/output folders based on system type
  if(path == "relative"){
    if(.Platform$OS.type == "unix") {
      path_db <- "~/Dropbox"
      path_raw <- "~/Dropbox/eutops/data/raw_data/methylight/"
      path_processed <- "~/Dropbox/eutops/data/processed_data/methylight/"
      
    } else {
      
      path_db <- path_expand("~/Dropbox")
      path_raw <- paste0(path_db, "/eutops/data/raw_data/methylight/")
      path_processed <- paste0(path_db, "/eutops/data/processed_data/methylight/")
    }
  
    # path to input folder
    folder <- paste0(path_raw, folder)
    
    # path to output folder
    output <- paste0(path_processed, output)
    
  }
  
  # Get scripts
  if(.Platform$OS.type == "unix"){
    source("~/Dropbox/eutops/scripts/methylight/pmr_calc/src/calculate_pmr.R")
    autotest_folder <- "~/Dropbox/eutops/scripts/methylight/pmr_calc/src/autotest-plate/"
    autotest_output <- "~/Dropbox/eutops/scripts/methylight/pmr_calc/src/tmp/"
    autotest.dat <- readRDS("~/Dropbox/eutops/scripts/methylight/pmr_calc/src/autotest-data.Rds")
    
  } else {
    path_db <- path_expand("~/Dropbox")
    source(paste0(path_db, "/eutops/scripts/methylight/pmr_calc/src/calculate_pmr.R"))
    autotest_folder <- paste0(path_db, "/eutops/scripts/methylight/pmr_calc/src/autotest-plate/")
    autotest_output <- paste0(path_db, "/eutops/scripts/methylight/pmr_calc/src/tmp/")
    autotest.dat <- readRDS(paste0(path_db, "/eutops/scripts/methylight/pmr_calc/src/autotest-data.Rds"))
  }
  
  
  cat("Data run beginning ...\n")
  
  results <- calculate_pmr(folder, output, experimentname,
                           write.results,
                           type,
                           calib_fixed, 
                           fix_intercept,
                           fix_slope,
                           threshold_COL2A1, path)
  
  cat("\nData run done ==>\nAUTOTEST beginning...\n")
  
  suppressWarnings(invisible(capture.output(autotest <- calculate_pmr(autotest_folder, autotest_output, "autotest",
                            write.results = F,
                            type = "gblock", calib_fixed = F, path = "fixed",
                            is.autotest = T) %>%
    filter(! Sample %in% c("H2O", "STD_1", "STD_2", "STD_3", "STD_4", "gBLOCK")))))
  
  if(!identical(autotest$Sample, autotest.dat$Sample)){
    
    if(sum(autotest$Sample %in% autotest.dat$Sample) != nrow(autotest.dat)){
      msg <- "WARNING: new script does not return all autotest samples and may remove certain samples from the output.\n"
      cat(msg)
      write(msg, file = paste0(output, experimentname, "_log.txt"), append = T)
      
      intersect <- intersect(autotest.dat$Sample, autotest$Sample)
      autotest.dat <- autotest.dat[match(intersect, autotest.dat$Sample),]
      autotest <- autotest[match(intersect, autotest$Sample),]
    }
    
    ind <- match(autotest.dat$Sample, autotest$Sample)
    autotest <- autotest[ind,]
    
  }
  
  autotest <- autotest %>%
    dplyr::select(colnames(autotest.dat))
  autotest$PMR_EC <- autotest$ZSCAN12 + autotest$GYPC1 + autotest$GYPC2
  autotest.dat$PMR_EC <- autotest.dat$ZSCAN12 + autotest.dat$GYPC1 + autotest.dat$GYPC2
  
  autotest1 <- autotest %>%
    mutate(set = "new") %>%
    pivot_longer(cols = -c("Sample", "set"),
                 names_to = "target",
                 values_to = "pmr")
  autotest.dat1 <- autotest.dat %>%
    mutate(set = "original") %>%
    pivot_longer(cols = -c("Sample", "set"),
                 names_to = "target",
                 values_to = "pmr")
  tmp <- rbind(autotest1, autotest.dat1) %>%
    pivot_wider(id_cols = c("Sample", "target"),
                names_from = "set", 
                values_from = "pmr")
  
  testplot <- tmp %>%
    ggplot(aes(x = log(original+0.000001),
               y = log(new+0.000001))) +
    geom_smooth(method = "lm",
                formula = y ~ x,
                size = 0.5) +
    geom_point(size = 0.5) +
    facet_wrap(~target) +
    stat_cor(method = "pearson") +
    theme_gray() +
    theme(panel.grid = element_blank()) +
    xlab("log(original values)") +
    ylab("log(current values)")
  
  ggsave(testplot,
         file = paste0(output, "autotest", ".png"),
         width = 6,
         height = 6) 
  
  cat("Autotest result: ",
      ifelse(isTRUE(identical(autotest, autotest.dat)), "PASS", "FAIL"))
  
  write(paste0("\n\n------------ AUTOTEST REPORT ------------\nAll expected samples in output: ",
               ifelse(isTRUE(sum(autotest$Sample %in% autotest.dat$Sample) == nrow(autotest.dat)), "Yes", "NO (!)"), "\nValues identical for all targets: ",
               ifelse(isTRUE(identical(autotest, autotest.dat)), "Yes", "NO (!)"),
               "\nAutotest result: ",
               ifelse(isTRUE(identical(autotest, autotest.dat)), "PASS", "FAIL"),
               "\n------------ RUN DONE ------------"),
        file = paste0(output, experimentname, "_log.txt"), append = T)
  
  
  ## Tidy up
  # Remove autotest plot
  f <- list.files(autotest_output, full.names = T)
  invisible(file.remove(f))

  return(results)
  
  }
  
  
  
  


