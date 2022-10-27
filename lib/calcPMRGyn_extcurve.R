# PMR calculation main work function
# contact: charlottevavourakis@gmail.com
# under development, Last update October 26, 2022
# this version has the correct way of calculating coefficients for multi-point calibration
# this version does not calculate the mean, until after threshold_COL2A1 checking is done
# Ct Sd is given again
# everything is calculated from the results file
# tested on the Barcelona data that the calibration correction does not make a difference in the final outcome of the WID-qEC (besides giving slightly different PMRs) compared to original pipeline published
# tested on one plate of the Hall dataset, that the fixed intercept/slope option works the same

# Option for using fixed intercept or slope, or external calibration curve, will only be functional in the bare R-pipeline, not in the Shiny apps


calcPMRGynExt <- function(samplesdata,gblockdata,intercept, slope,threshold_COL2A1, threshold_targets){
 
  require(dplyr)
  require(tidyverse)
  
  # initialize list so multiple objects can be returned from function ----
  listResults = list()
  
  ####---- Calculate gblocks separately ----####
  
  data <- gblockdata %>%
    mutate(rep = case_when(grepl("A|C|E|G|I|K|M|O", Well.Position) == TRUE ~ 1,
                           grepl("B|D|F|H|J|L|N|P", Well.Position) == TRUE ~ 2))
  
  # collect target and sample overviews
  targets <- unique(data$Target)
  targets2 <- targets[!targets %in% "COL2A1"]
  samples <- unique(data$Sample)
  
  ####---- Set suspicious amplifications to NA----####
  data <- data %>% 
    mutate(ct = as.numeric(ifelse(Cq == "Undetermined", NA, Cq.Mean))) %>% 
    mutate(ct = as.numeric(ifelse(Cq > threshold_targets, NA, Cq)))
  
  # If one rep is undetermined, set both to NA (undetermined)
  samples <- unique(data$Sample)
  targets <- unique(data$Target)
  reprocess <- ""
  
  for (i in 1:length(samples)){
    for(j in 1:length(targets)){
      tmp <- data %>%
        filter(Sample == samples[i] & Target == targets[j]) %>%
        pull(ct)
      
      und <- sum(is.na(tmp))
      ind1 <- data$Sample==samples[i] & data$Target==targets[j]
      ct <- ifelse(length(ind1)==1 && is.na(tmp), 34, ifelse((length(ind1)==1 && !is.na(tmp)), tmp, ifelse(und<2, na.omit(tmp), NA)))
      
      if(sum(ind1)==1 && is.na(tmp)){
        data[ind1,]$ct <- NA
      } else if(sum(ind1)==1 && !is.na(tmp)){
        data[ind1,]$ct <- ct
      } else if(und > 0){
        reprocess <- c(reprocess, data[ind1,]$Sample)
        data[ind1,]$ct <- rep(ct, 2)
      }
    }
  }
  
  reprocess <- unique(reprocess[-1])
  
  data <- data %>%
    pivot_wider(id_cols = Sample,
                names_from = Target,
                values_from = ct,
                values_fn = mean)
  
  # Calculate input amounts
  data <- data %>%
    mutate(input_ZSCAN12 = 10^((ZSCAN12-intercept)/slope),
           input_GYPC1 = 10^((GYPC1-intercept)/slope),
           input_GYPC2 = 10^((GYPC2-intercept)/slope),
           input_DPP6 = 10^((DPP6-intercept)/slope),
           input_RALYL = 10^((RALYL-intercept)/slope),
           input_GSX1 = 10^((GSX1-intercept)/slope),
           input_COL2A1 = 10^((COL2A1-intercept)/slope),
           input_EpC = 10^((EpC-intercept)/slope),
           input_ImC = 10^((ImC-intercept)/slope)) %>%
    mutate(ref_ZSCAN12 = input_ZSCAN12/input_COL2A1,
           ref_GYPC1 = input_GYPC1/input_COL2A1,
           ref_GYPC2 = input_GYPC2/input_COL2A1,
           ref_DPP6 = input_DPP6/input_COL2A1,
           ref_RALYL = input_RALYL/input_COL2A1,
           ref_GSX1 = input_GSX1/input_COL2A1,
           ref_COL2A1 = input_COL2A1/input_COL2A1,
           ref_EpC = input_EpC/input_COL2A1,
           ref_ImC = input_ImC/input_COL2A1)
  
  
  # Separate for gblocks
  gblock <- data %>%
    filter(Sample == "gBLOCK") 
  
  ####---- Calculate samples ----####
  ####---- Set suspicious amplifications to NA----####
  
  data = samplesdata %>% mutate(ct = as.numeric(ifelse(Cq == "Undetermined", NA, Cq)))
  #data = data %>% mutate(ct= as.numeric(ifelse(Amp.Status != "Amp", NA, Cq))) #this option was voted out
  data = data %>% mutate(ct = as.numeric(ifelse(Cq > threshold_targets, NA, Cq)))
  
  # collect target and sample overviews
  targets <- unique(data$Target)
  targets2 <- targets[!targets %in% "COL2A1"]
  samples <- unique(data$Sample)
  
  ####---- Check COL2A1 ----####
  # For each rep evaluate COL2A1 and add a column PASS/FAIL for filtering data
  df = NULL
  for (s in samples) {
    subset1 = data %>%
      filter(Sample==s,
             rep==1,
             Target=="COL2A1")
    subset1 = subset1 %>%
      mutate(COL2A1_check =
               case_when( 
                 is.na(ct) | ct > threshold_COL2A1 ~ "FAIL", TRUE ~ "PASS")) %>%
      dplyr::select(Sample,rep,COL2A1_check)
    df = rbind(df,subset1)
    
    subset2 = data %>%
      filter(Sample==s,
             rep==2,
             Target=="COL2A1")
    subset2 = subset2 %>%
      mutate(COL2A1_check = 
               case_when( 
                 is.na(ct) | ct > threshold_COL2A1 ~ "FAIL", TRUE ~ "PASS")) %>%
      dplyr::select(Sample,rep,COL2A1_check)
    df = rbind(df,subset2)
  }
  
  data = left_join(data, df, by=c("Sample","rep"))
  
  # When both reps of a sample are undetermined or above the set threshold for COL2A1,
  # DNA input for this sample was too low, and no PMR could be calculated.
  # If one out of two reps is undetermined or above the set threshold,
  # sample should to be reprocessed (COL2A1) or reprocessing is recommended. 
  # In these cases, one Ct value is taken forward. 
  low_input_fail = NULL # negative control should always endup here
  reprocess_needed = NULL
  reprocess = NULL
  
  for (s in samples){
    for(t in targets){
      tmp = data %>%
        filter(Sample == s & Target == t) %>%
        dplyr::select(Sample, rep, ct, COL2A1_check)
      
      if(t == "COL2A1" & all(tmp$COL2A1_check=="FAIL")){
        # if target is COL2A1 and both reps fail, sample should be listed under low_input_fail
        x = data.frame(sample = tmp$Sample, COL2A1_Ct = tmp$ct)
        low_input_fail = rbind(low_input_fail, x)
      }
      
      if(any(is.na(tmp$ct)) && any(!is.na(tmp$ct))){
        # Carry one Ct value forward
        und = sum(is.na(tmp$ct))
        ind1 = na.omit(data$Sample == s & data$Target == t)
        ct = ifelse(und < 2, na.omit(tmp$ct), NA)
        data[ind1,]$ct <- ct
        
        if(t == "COL2A1"){
          # reprocessing required
          x = data.frame(sample = tmp$Sample[is.na(tmp$ct)],
                         rep = tmp$rep[is.na(tmp$ct)],
                         target = t,
                         `COL2A1_check (this replicate)` = tmp$COL2A1_check[is.na(tmp$ct)])
          reprocess_needed = rbind(reprocess_needed,x)
        } else {
          # reprocessing recommended
          x = data.frame(sample = tmp$Sample[is.na(tmp$ct)],
                         rep = tmp$rep[is.na(tmp$ct)],
                         target = t,
                         `COL2A1_check (this replicate)` = tmp$COL2A1_check[is.na(tmp$ct)])
          reprocess = rbind(reprocess,x)
        }
      }
    }
  }
  
  data1 = data %>%
    pivot_wider(id_cols = Sample,
                names_from = Target,
                values_from = ct,
                values_fn = mean)
  
  data2 = data %>%
    pivot_wider(id_cols = Sample,
                names_from = Target,
                values_from = ct,
                values_fn = sd)
  
  # Calculate PMR ----
  
  # Calculate input amounts (using calibration curve COL2A1)
  for (t in targets){
    varname = paste0("input_",t)
    data1[varname] = with(data1, 10^((data1[t]-intercept)/slope))
  }
  
  # Normalization: Divide calculated input amount by input amount calculated for COL2A1
  for (t in targets){
    varname = paste0("ref_",t)
    calcname = paste0("input_",t)
    data1[varname] = data1[calcname]/data1["input_COL2A1"]
  }
  
  # Fully methylated DNA (SSS1 or gBlock) - only gBlock implemented
  gblock = data1 %>%
    filter(Sample %in% c("posCo","PosCo","gBlock", "gBLOCK", "gBlock (+)")) 
  
  # For each ref gene divide normalized input sample by normalized input fully methylated DNA *100
  results = data1
  
  for (t in targets2){
    varname = t
    calcname = paste0("ref_",t)
    div_block = c(gblock[calcname])
    results[varname] = (data1[calcname]/div_block)*100
  }
  
  results = as.data.frame(results)
  results = results %>%
    select(Sample, all_of(targets2))
  results[is.na(results)] = 0 # set all NAs to 0
  
  # Calculate concentration input DNA based on COL2A mean Ct
  conc_input = data1 %>% 
    as.data.frame() %>%
    mutate(input_conc = (data1$COL2A1 - intercept)/slope) %>%
    select(Sample, input_conc) %>% droplevels
  
  # Calculate WID-qEC and or WID-CIN
  targets_qEC = c("GYPC1","GYPC2","ZSCAN12")
  if(all(targets_qEC %in% targets)){
    results = results %>%
      mutate(WIDqEC = round((GYPC1 + GYPC2 + ZSCAN12), digits=3)) # calculate sumPMR for WID-qEC
  }
  
  targets_qCIN = c("RALYL", "DPP6", "GSX1")
  if(all(targets_qCIN %in% targets)){
    results = results %>%
      mutate(WIDqCIN = round((RALYL + DPP6 + GSX1), digits=3)) # calculate sumPMR for WID-qCIN
  }
  
  # samples for which both reps COL2A1 failed, PMR should not be 0, but should be NA
  if(!is_empty(low_input_fail)){
    samples_failed = low_input_fail %>% filter(!sample %in% c("NTC","NTC_H2O","H2O")) %>% unique() %>% pull(sample)
    if(!is_empty(samples_failed)){
      for (i in 2:(length(targets2)+3)){
        for (j in 1:length(samples_failed)){
          results[results$Sample==samples_failed[j],i] = NA
        }
      }
    }
  }
  
  # Save results to list ---
  data1 = data1 %>%
    select(Sample, COL2A1, all_of(targets2)) %>%
    as.data.frame()
  
  data2 = data2 %>%
    select(Sample, COL2A1, all_of(targets2)) %>%
    as.data.frame()
  
  listResults[[1]] = results #PMR + optional WIDqEC + WIDqEC outcome EUTOPS + WIDqEC interpret TP +WIDqCIN
  listResults[[2]] = data1 #mean Cq
  listResults[[3]] = data2 #stdev Cq
  listResults[[4]] = conc_input #log(copy number/5uL)
  if(is_empty(low_input_fail)){
    listResults[[5]] = as.data.frame("you have some amplification in the NTC_H2O")
  } else{
    listResults[[5]] = low_input_fail #samples for which COL2A1 failed in both reps, should have negative controls
  }
  
  if(!is_empty(reprocess_needed)){ 
    listResults[[6]] = reprocess_needed
  } else {
    listResults[[6]] = as.data.frame("none")
  }
  
  if(!is_empty(reprocess)){ # samples for which for only one of two reps target amplified
    listResults[[7]] = reprocess
  } else {
    listResults[[7]] = as.data.frame("none")
  }
  
  return(listResults)
}
