# wraps COL2A1 calibration, PMR calculation and testWIDqEC main work functions, based on user input
# Option for using fixed intercept or slope, or external calibration curve, will only be functional in the bare R-pipeline, not in the Shiny apps

calculate_pmr <- function(data,
                          threshold_COL2A1,
                          threshold_targets,
                          external_curve,
                          fix_intercept,
                          fix_slope,
                          calib_fixed,
                          inconclusive, 
                          inconsistent, 
                          version){
  
  ####---- Prepare ----####
  
  require(dplyr)
  require(tidyverse)
  
  # main work functions
  source("../../lib/calibGblock.R")
  source("../../lib/calcPMRGyn.R")
  source("../../lib/testWIDqEC.R")
  source("../../lib/plotPlateCTmeanSD.R")
  
  # check that threshold targets >= threshold COL2A1
  if(threshold_COL2A1>threshold_targets){
    stop("Error: CT-threshold set for COL2A1 should be lower or equal to that for other targets")
  }
  
  # collect target and sample overviews
  targets <- unique(data$Target)
  samples <- unique(data$Sample)
  samplesNoctrl <- setdiff(samples,c("STD1", "STD2", "STD3","STD4","STD_1", "STD_2", "STD_3", "STD_4","Std 1", "Std 2", "Std 3", "Std 4","Std_1", "Std_2", "Std_3", "Std_4","Std. 1", "Std. 2", "Std. 3", "Std. 4","NTC_H2O","NTC", "H2O")) #but include gBlocks C("posCo","PosCo","gBlock", "gBLOCK", "gBlock (+)")
  
  # initialize list so multiple objects can be returned from function 
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
  
  ####----- Plot mean+SD CT COL2A1 ----####
  
  list_out[[10]] <- plotCTPlate(samplesNoctrl,"COL2A1",calcPMR[[2]],calcPMR[[3]])
  
  ####---- Make the final summary commercial WIDqEC test ----####
  
  targets_qEC <- c("GYPC1","GYPC2","ZSCAN12")
  if(all(targets_qEC %in% targets)){
    testEC <- testWIDqEC(calcPMR[[1]],samples,calcPMR[[5]],calcPMR[[6]],calcPMR[[7]],calcPMR[[8]])
    list_out[[11]] <- testEC[[1]] #information in final csv file
    list_out[[2]] <- testEC[[2]] #update final batch results
  }
  
  ###---- Create log file ----####
  L4 <- ifelse(!is.null(external_curve),"External standard curve used",ifelse(calib_fixed==TRUE,"Fixed regression parameters used","on-plate-calibration"))
  L5 <- ifelse(inconclusive==FALSE, "Inconclusive amplification not disregarded","Inconclusive amplifications disregarded")
  L6 <- ifelse(inconsistent==FALSE,"All other amplifications carried forward", "Inconsistent targets rejected")
  loginfo <- c(version,
           paste0("CT threshold COL2A1 = ", threshold_COL2A1),
           paste0("CT threshold all other targets = ", threshold_targets),
           L4,
           L5,
           L6)
  list_out[[11]] <- as.data.frame(loginfo)
  
  ####---- Return or save results ----####
  
  return(list_out)
  
}
