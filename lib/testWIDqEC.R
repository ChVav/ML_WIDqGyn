## Summarize results WIDqEC test
# interpretation based on set thresholds
# sample QC

testWIDqEC <- function(results,samples,low_input_fail,reprocess_needed,reprocess){
  
  require(dplyr)
  require(tidyverse)
  
  # WID-qEC thresholds, set https://doi.org/10.1002/ijc.34275
  qEC_threshold1 <- 0.03
  qEC_threshold2 <- 0.63 #this threshold is not used by Tyrolpath
  
  # make the final summary diagnosis commercial WIDqEC test
  results = results %>% 
    mutate(WIDqEC_test = case_when(results$WIDqEC < qEC_threshold1 ~ "low risk EC/CIN",
                                     results$WIDqEC >= qEC_threshold1 & results$WIDqEC <= qEC_threshold2 ~ "high risk EC/CIN",
                                     results$WIDqEC > qEC_threshold2 ~ "very high risk EC/CIN")) %>%
    mutate(WIDqEC_interpret = case_when(results$WIDqEC >= qEC_threshold1 ~ "Positiv",
                                          TRUE ~ "Negativ"))
  final <- results %>%
    filter(!Sample %in% c("Std_1", "Std_2", "Std_3", "Std_4",
                            "Std. 1", "Std. 2", "Std. 3", "Std. 4",
                            "Std 1", "Std 2", "Std 3", "Std 4",
                            "STD1", "STD2", "STD3","STD4",
                            "STD_1","STD_2","STD_3","STD_4",
                            "posCo","PosCo","gBlock", "gBLOCK", "gBlock (+)",
                            "NTC_H2O","NTC", "H2O")) %>%
    select(Sample, WIDqEC, WIDqEC_interpret)
  final$WIDqEC <- format(final$WIDqEC, nsmall=3, scientific=FALSE)
  
  samples <- samples[!samples %in% c("Std_1", "Std_2", "Std_3", "Std_4",
                                       "Std. 1", "Std. 2", "Std. 3", "Std. 4",
                                       "Std 1", "Std 2", "Std 3", "Std 4",
                                       "STD1", "STD2", "STD3","STD4",
                                       "STD_1","STD_2","STD_3","STD_4",
                                       "posCo","PosCo","gBlock", "gBLOCK", "gBlock (+)",
                                       "NTC_H2O","NTC", "H2O")]
    
  # add QC
  QC <- data.frame(Sample = samples)
  
  if(!is_empty(low_input_fail)){
    QC <- QC %>%
      mutate(QC=case_when(
        Sample %in% low_input_fail$sample ~ "Insufficient DNA", #COL2A1 did not amplify in any of the reps 
        Sample %in% reprocess_needed[,1] ~ "Reprocessing recommended, insufficient DNA in one rep", # samples for which for only one of two reps COL2A1 failed
        Sample %in% reprocess[,1] ~ "Some targets only amplified in one of the reps",
        TRUE ~ "PASS" #all ok
      ))} else{
        QC <- QC %>%
          mutate(QC=case_when(
            Sample %in% reprocess_needed[,1] ~ "Reprocessing recommended, insufficient DNA in one rep", # samples for which for only one of two reps COL2A1 failed
            Sample %in% reprocess[,1] ~ "Some targets only amplified in one of the reps",
            TRUE ~ "PASS" #all ok
          ))
      }
  
  final <- full_join(final,QC)
  
  list <- list(final, results) #also results are updated
  
  return(list)
}
