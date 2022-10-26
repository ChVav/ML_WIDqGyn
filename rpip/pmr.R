# wraps COL2A1 calibration, PMR calculation and testWIDqEC main work functions, based on user input
# contact: charlottevavourakis@gmail.com
# R pipeline under development for automating analysis GYN test, Last update October 26, 2022
#
# Description:
#   This version works for WID-qGYN, with minimally COL2A1, EpC and ImC targets. Targets to exclude
#   are set in the previous part of the pipeline (sheet generation).
#   This is a wrapper function. It calls the most recent PMR script and calculates PMRs 
#
# Variables:
#   - calib_fixed: Instead of using the on-plate calibration curve, fixed values for intercept and slope
#     can be used.
#   - external_curve: Instead of using the on-plate calibration curve, an external calibration curve (run on a different plate may be used)
#   - path: relative (default: location in raw data folder) or fixed (input from other locations)
#   - autotest: perform autotest (default=FALSE)
#
# Still to be implemented:
# works only for gblock standards, sss1 may be implemented in the future
# autotest for the most up to date script to see whether it replicates previous
#   version's values on a barcelona plate.

pmr <- function(folder=NULL, output=NULL, experimentname,
                write.results = TRUE,
                threshold_COL2A1=30,
                threshold_targets=38,
                external_curve=NULL,
                calib_fixed = FALSE,
                fix_intercept = 36.9,
                fix_slope = -3.4,
                path = "relative",
                autotest= FALSE,
                path_lib="/eutops/scripts/methylight/ML_WIDqGYn/lib/"){
  
  ####---- Prepare ----####
  
  # load libraries
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(openxlsx))
  suppressPackageStartupMessages(library(fs))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggpubr))
  
  # specify paths input/output folders based on system type
  if(.Platform$OS.type == "unix") {
    path_db <- "~/Dropbox"
  } else {
    path_db <- path_expand("~/Dropbox")
  }
  
  if(path == "relative"){
    path_raw <- paste0(path_db, "/eutops/data/raw_data/methylight/")
    path_processed <- paste0(path_db, "/eutops/data/processed_data/methylight/")
    folder <- paste0(path_raw, folder)  # path to input folder
    output <- paste0(path_processed, output) # path to output folder
  }
  
  # source main work functions
  if(path_lib == "/eutops/scripts/methylight/ML_WIDqGYn/lib/"){
    path_lib <- paste0(path_db,path_lib)
  }
  
  source(paste0(path_lib,"calibGblock.R"))
  source(paste0(path_lib,"calcPMRGyn.R"))
  source(paste0(path_lib,"testWIDqEC.R"))
  
  # check that threshold targets >= threshold COL2A1
  if(threshold_COL2A1>threshold_targets){
    stop("Error: CT-threshold set for COL2A1 should be lower or equal to that for other targets")
  }
  
  # find files: if required input files not present, stop; if files in output directory already exists, throw warning (preventing overwriting)
  files <- list.files(folder)
  #curve <- grep("Standard Curve", files)
  results <- grep("Results", files)
  
  #if(length(curve) == 0){
  #  stop("Error: file with standard curves not found")
  #}
  
  if(length(results) == 0){
    stop("Error: results file not found in provided folder")
  }
  
  if(dir.exists(output)){
    cat(paste0("Output folder ", output, " already exists.\n"))
  } else {
    dir.create(output, recursive = TRUE)
    cat(paste0("Output folder ", output, " created.\n"))
  }
  
  output_folder <- list.files(output)
  if(paste(experimentname, ".xlsx", sep = "") %in% output_folder == TRUE && write.results == TRUE){
    stop("File with the same name already exists in the output folder. Please choose a unique name.")
  }
  
  # Startup messages
  cat("PMR calculation starting...\n")
  
  msg2 <- paste0("[pmr] Settings:\n [pmr] folder = ", folder, "\n [pmr] output = ", output, "\n [pmr] experiment name = ", experimentname, "\n [pmr] threshold_COL2A1 = ", threshold_COL2A1, "\n [pmr] threshold_targets = ", threshold_targets, "\n [pmr] external_curve = ", external_curve, ifelse(external_curve!=NULL, paste0(" (CAUTION: using external calibration curve)"), ""), "\n [pmr] calib_fixed = ", calib_fixed, ifelse(isTRUE(calib_fixed), paste0(" (CAUTION: using fixed intercept/slope values)"), ""), "\n [pmr] path = ", path, sep = "","\n")

  cat(msg2)
  
  # Log message (including warning)
  write(paste0("[pmr] ", date()), file = paste0(output, experimentname, "_log.txt"), append = F)
  write(paste0("[pmr] Platform: ", .Platform$OS.type), file = paste0(output, experimentname, "_log.txt"), append = T)
  write(msg2, file = paste0(output, experimentname, "_log.txt"), append = T)
  
  if(calib_fixed == TRUE){
    write("\nWarning: instead of using the on-plate calibration curve, you are using FIXED values for the intercept and slope. Use this option with caution.\n",
          file = paste0(output, experimentname, "_log.txt"), append = T)
  }
  
  # read in data
  data <- read.table(file = paste(folder, files[results], sep = ""),
                     sep = ",",
                     header = TRUE)
  
  # collect target and sample overviews
  targets <- unique(data$Target)
  samples <- unique(data$Sample)
  
  # initialize list so multiple objects can be returned from main work functions
  list_out <- list()
  
  ####---- Calibrate COL2A1 ----####
  
  calibration <- calibGblock(data)
  
  # add plot to list
  list_out[[1]] <- calibration[[1]]
  
  # If user choose fixed values, overwrite the coefficients from standard curve
  if(calib_fixed==TRUE){
    intercept <- fix_intercept
    slope <- fix_slope
  } else{
    intercept <- calibration[[2]]
    slope <- calibration[[3]]
  }
  
  ####---- Reformat data ~reps knowing plate layout  ----####
  
  data <- data %>%
    mutate(rep = case_when(grepl("A|C|E|G|I|K|M|O", Well.Position) == TRUE ~ 1,
                           grepl("B|D|F|H|J|L|N|P", Well.Position) == TRUE ~ 2))
  
  ####---- Calculate PMR ----####
  # based on default or user set thresholds
  # minimum two reps needed, not sure whether works for more
  
  calcPMR <- calcPMRGyn(data,targets,samples,intercept, slope,threshold_COL2A1, threshold_targets)
  list_out <- append(list_out,calcPMR)
  
  ####---- Make the final summary commercial WIDqEC test ----####
  
  targets_qEC <- c("GYPC1","GYPC2","ZSCAN12")
  if(all(targets_qEC %in% targets)){
    testEC <- testWIDqEC(calcPMR[[1]],samples,calcPMR[[5]],calcPMR[[6]],calcPMR[[7]])
    list_out[[9]] <- testEC[[1]] #information in final csv file
    list_out[[2]] <- testEC[[2]] #update final batch results
  }

  ####---- Autotest ---- ####
  
  ## note this part is not working yet for current version
  if(autotest==TRUE){
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
  }
  
  ####---- Return and/or write results ----####
  
  if(write.results == TRUE){
    
    ggsave(list_out[[1]],
           file = paste0(output,experimentname, ".png"),
           width = 5,
           height = 4) 
    cat("COL2A1 calibration curve saved under", paste0(output,experimentname, ".png\n"))
    
    wb <- createWorkbook()
    addWorksheet(wb, "Results PMR")
    writeDataTable(wb = wb, sheet = 1, x = list_out[[2]], rowNames=FALSE)
    addWorksheet(wb, "Ct means")
    writeDataTable(wb = wb, sheet = 2, x = list_out[[3]], rowNames=FALSE)
    addWorksheet(wb, "Ct sd")
    writeDataTable(wb = wb, sheet = 3, x = list_out[[4]], rowNames=FALSE)
    addWorksheet(wb, "Input BC-DNA conc")
    writeDataTable(wb=wb, sheet = 4, x= list_out[[5]], rowNames=FALSE)
    addWorksheet(wb, "low DNA input") #samples for which COL2A1 failed in both reps, should have negative controls
    if(is_empty(list_out[[6]])){
      low_input_fail <- as.data.frame("you have some amplification in the NTC_H2O")
      writeDataTable(wb=wb, sheet = 5, x= list_out[[6]], rowNames=FALSE)
    } else{
      writeDataTable(wb=wb, sheet = 5, x= list_out[[6]], rowNames=FALSE) 
    }
    addWorksheet(wb, "Reprocessing needed")# samples for which for only one of two reps COL2A1 failed
    if(!is_empty(list_out[[7]])){
      writeDataTable(wb = wb, sheet = 6, x = list_out[[7]], rowNames=FALSE)
    } else{
      reprocess_needed <- as.data.frame("none")
      writeDataTable(wb = wb, sheet = 6, x = list_out[[7]], rowNames=FALSE)
    }
    addWorksheet(wb, "Reprocessing recommended")# samples for which for only one of two reps target amplified
    if(!is_empty(list_out[[8]])){
      writeDataTable(wb = wb, sheet = 7, x = list_out[[8]], rowNames=FALSE)
    } else{
      reprocess <- as.data.frame("none")
      writeDataTable(wb = wb, sheet = 7, x = list_out[[8]], rowNames=FALSE)
    }
    
    saveWorkbook(wb, paste0(output, experimentname, ".xlsx"))
    
    cat("Output file created.")
  }

  return(list_out)
  
  }
  
  
  
  

