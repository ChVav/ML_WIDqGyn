# PMR calculation main work function

calculate_pmr <- function(folder, output, experimentname,
                          write.results = TRUE,
                          type = "gblock",
                          calib_fixed = FALSE,
                          fix_intercept,
                          fix_slope,
                          threshold_COL2A1,
                          path = "relative", 
                          is.autotest = FALSE){
  
  # find files: if required input files not present, stop; if files in output directory already exists, throw warning (preventing overwriting)
  files <- list.files(folder)
  curve <- grep("Standard Curve", files)
  results <- grep("Results", files)
  
  if(length(curve) == 0){
    stop("Error: file with standard curves not found")
  }
  
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
  if(isFALSE(is.autotest)){
    msg1 <- ("[pmr_calc] PMR calculation starting...\n")
  } else {
    msg1 <- ("[pmr_calc] PMR calculation for automatic test run starting...\n")
  }
  msg2 <- paste0("[pmr_calc] Settings:\n[pmr_calc]    folder = ", folder, "\n[pmr_calc]    output = ", output, "\n[pmr_calc]    experiment name = ", experimentname, "\n[pmr_calc]    type = ", type, "\n[pmr_calc]    calib_fixed = ", calib_fixed, ifelse(isTRUE(calib_fixed), paste0(" (CAUTION: using fixed intercept/slope values)"), ""), "\n[pmr_calc]    threshold_COL2A1 = ", threshold_COL2A1, "\n[pmr_calc]    path = ", path, sep = "")
  
  cat(msg1)
  cat(msg2)
  
  # Log message (including warning)
  write(paste0("[pmr_calc] ", date()), file = paste0(output, experimentname, "_log.txt"), append = F)
  write(paste0("[pmr_calc] Platform: ", .Platform$OS.type), file = paste0(output, experimentname, "_log.txt"), append = T)
  write(msg2, file = paste0(output, experimentname, "_log.txt"), append = T)

  
  if(calib_fixed == TRUE){
    write("\nWarning: instead of using the on-plate calibration curve, you are using FIXED values for the intercept and slope. Use this option with caution.\n",
                file = paste0(output, experimentname, "_log.txt"), append = T)
  }
  
  
  # Generate standard curve from COL2A1 values
  curve <- read.table(file = paste(folder, files[curve], sep = ""),
                      sep = ",",
                      header = TRUE)
  
  curve <- suppressWarnings(curve %>%
                              mutate(ct = as.numeric(Cq.Mean)))
  
  if(type == "gblock"){
    plot <- curve %>%
      filter(Target == "COL2A1" & ((Sample %in% c("STD1", "STD2", "STD3","STD4") | (Sample %in% c("STD_1", "STD_2", "STD_3", "STD_4")) | (Sample %in% c("Std 1", "Std 2", "Std 3", "Std 4")) | (Sample %in% c("Std_1", "Std_2", "Std_3", "Std_4")) | (Sample %in% c("Std. 1", "Std. 2", "Std. 3", "Std. 4"))))) %>%
      mutate(x = case_when(Sample %in% c("STD1", "STD_1", "Std 1", "Std_1", "Std. 1") ~ 6,
                           Sample %in% c("STD2", "STD_2", "Std 2", "Std_2", "Std. 2") ~ 5,
                           Sample %in% c("STD3", "STD_3", "Std 3", "Std_3", "Std. 3") ~ 4,
                           Sample %in% c("STD4", "STD_4", "Std 4", "Std_4", "Std. 4") ~ 3)) %>%
      mutate(ct = ifelse(Cq == "Undetermined", NA, as.numeric(Cq)))
  } else {
    plot <- curve %>%
      filter(Target == "COL2A1" & ((Sample %in% c("STD1", "STD2", "STD3","STD4") | (Sample %in% c("STD_1", "STD_2", "STD_3", "STD_4")) | (Sample %in% c("Std 1", "Std 2", "Std 3", "Std 4")) | (Sample %in% c("Std_1", "Std_2", "Std_3", "Std_4")) | (Sample %in% c("Std. 1", "Std. 2", "Std. 3", "Std. 4"))))) %>%
      mutate(x = case_when(Sample %in% c("STD1", "STD_1", "Std 1", "Std_1", "Std. 1") ~ 250,
                           Sample %in% c("STD2", "STD_2", "Std 2", "Std_2", "Std. 2") ~ 62.5,
                           Sample %in% c("STD3", "STD_3", "Std 3", "Std_3", "Std. 3") ~ 15.625,
                           Sample %in% c("STD4", "STD_4", "Std 4", "Std_4", "Std. 4") ~ 3.9)) %>%
      mutate(ct = ifelse(Cq == "Undetermined", NA, as.numeric(Cq))) 
  }
  
  fit <- lm(ct ~ x, data = plot)
  anno <- paste0("R2 = ", signif(summary(fit)$r.squared, 3),
                 "\ny = ", round(fit$coefficients[2],2) , "x + ", round(fit$coefficients[1],2))
  y = max(plot$ct) - 0.05*max(plot$ct)
  
  # Warning low R.squared
  R2 <- summary(fit)$r.squared
  if (R2<0.99 && calib_fixed==FALSE){
    write(paste0("Calibration curve R2 = ", signif(R2, 4), "\nWARNING: Your COL2A1 calibration curve is below the ideal of 0.99. You may consider rerunning the PMR calculation using fixed values for the calibration intercept and slope by setting the option calib_fixed = TRUE.", sep = ""), file = paste0(output, experimentname, "_log.txt"), append = T)
          
  } else if (calib_fixed==TRUE){
    write(paste0("Calibration curve R2 = ", signif(R2, 4), "\nWARNING: Your on-plate calibration curve was not used for any calculations as you chose to use a fixed threshold.\n", sep = ""), file = paste0(output, experimentname, "_log.txt"), append = T)
  } else {
    write(paste0("Calibration curve R2 = ", signif(R2, 4)), file = paste0(output, experimentname, "_log.txt"), append = T)
  }
  
  # Plot calibration curve
  curveplot <- plot %>%
    ggplot(aes(x = x,
               y = ct)) +
    geom_smooth(method = "lm",
                formula = y ~ x,
                size = 0.5,
                se = FALSE,
                alpha = 0.3,
                colour = "gray40") +
    geom_point(aes(colour = Sample),
               size = 0.75) +
    theme(legend.title = element_blank(),
          legend.position = "top",
          legend.key = element_blank()) +
    xlab("Concentration") +
    ylab("Ct value")  +
    annotate("text",
             x = 5,
             y = 28,
             label = anno,
             hjust = 0,
             size = 3,
             colour = ifelse(R2 < 0.99, "red", "black"))
  
  if(isFALSE(is.autotest)){
  suppressMessages(print(curveplot))
  
  ggsave(curveplot,
         file = paste0(output,experimentname, ".png"),
         width = 5,
         height = 4) 
  write(paste0("COL2A1 calibration curve saved under ", output, experimentname, ".png\n"),
        file = paste0(output, experimentname, "_log.txt"), append = T)
  }
  
  rm(plot, curveplot)
  
  if(calib_fixed==FALSE){
    model <- curve %>%
      filter(Target == "COL2A1" & ((Sample %in% c("STD1", "STD2", "STD3","STD4") | (Sample %in% c("STD_1", "STD_2", "STD_3", "STD_4")) | (Sample %in% c("Std 1", "Std 2", "Std 3", "Std 4")) | (Sample %in% c("Std_1", "Std_2", "Std_3", "Std_4")) | (Sample %in% c("Std. 1", "Std. 2", "Std. 3", "Std. 4"))))) %>%
      group_by(Sample) %>%
      summarise(ct = mean(ct),
                quantity = log10(mean(Quantity))) %>%
      mutate(quantity = case_when(Sample %in% c("STD_1", "STD1", "Std 1", "Std_1", "Std. 1") ~ 6,
                                  Sample %in% c("STD_2", "STD2", "Std 2", "Std_2", "Std. 2") ~ 5,
                                  Sample %in% c("STD_3", "STD3", "Std 3", "Std_3", "Std. 3") ~ 4,
                                  Sample %in% c("STD_4", "STD4", "Std 4", "Std_4", "Std. 4") ~ 3))
    
    std <- lm(ct ~ quantity, data = model)
    intercept <- unname(std$coefficients[1])
    slope <- unname(std$coefficients[2])
  } else { 
    intercept <- fix_intercept
    slope <- fix_slope
  }
  
  # Get and format data
  data <- read.table(file = paste(folder, files[results], sep = ""),
                     sep = ",",
                     header = TRUE)
  
  data <- data %>%
    mutate(rep = case_when(grepl("A|C|E|G|I|K|M|O", Well.Position) == TRUE ~ 1,
                           grepl("B|D|F|H|J|L|N|P", Well.Position) == TRUE ~ 2)) %>%
    mutate(ct = as.numeric(ifelse(Cq == "Undetermined", NA, Cq.Mean)))
  
  targets <- unique(data$Target)
  targets2 <- targets[!targets %in% "COL2A1"]
  # standards <- c("STD_1", "STD1", "Std 1", "Std_1", "Std. 1","STD_2", "STD2", "Std 2", "Std_2", "Std. 2","STD_3", "STD3", "Std 3", "Std_3", "Std. 3","STD_4", "STD4", "Std 4", "Std_4", "Std. 4")
  # 
  # data <- data %>%
  #   filter(!Sample %in% standards) # retain only samples, negative (H2O and positive control)
  samples <- unique(data$Sample)
  
  # For each rep evaluate COL2A1 and add a column PASS/FAIL for filtering data
  # this does not make much sense when the mean is already taken, COL2A1 is not checked for the reps individually
  df <- NULL
  for (s in samples) {
    subset1 <- data %>%
      filter(Sample==s,
             rep==1,
             Target=="COL2A1")
    subset1 <- subset1 %>%
      mutate(COL2A1_check =
               case_when( 
      is.na(ct) | ct > threshold_COL2A1 ~ "FAIL", TRUE ~ "PASS")) %>%
      dplyr::select(Sample,rep,COL2A1_check)
    df <- rbind(df,subset1)
    
    subset2 <- data %>%
      filter(Sample==s,
             rep==2,
             Target=="COL2A1")
    subset2 <- subset2 %>%
      mutate(COL2A1_check = 
               case_when( 
      is.na(ct) | ct > threshold_COL2A1 ~ "FAIL", TRUE ~ "PASS")) %>%
      dplyr::select(Sample,rep,COL2A1_check)
    df <- rbind(df,subset2)
  }
  
  if(unique(df[df$Sample=="H2O",]$COL2A1_check) != "FAIL"){
    write("WARNING: COL2A1 amplified in the negative control.", file = paste0(output, experimentname, "_log.txt"), append = T)
  }
  
  data <- left_join(data, df, by=c("Sample","rep"))
  
  # When both reps of a sample are undetermined or above the set threshold for COL2A1,
  # DNA input for this sample was too low, and no PMR could be calculated.
  # If one out of two reps is undetermined or above the set threshold,
  # sample should to be reprocessed (COL2A1) or reprocessing is recommended. 
  # In these cases, one Ct value is taken forward. 
  low_input_fail <- NULL # negative control should always endup here
  reprocess_needed <- NULL
  reprocess <- NULL
  
  for (s in samples){
    for(t in targets){
      tmp <- data %>%
        filter(Sample == s & Target == t) %>%
        dplyr::select(Sample, rep, ct, COL2A1_check)
      
      if(t == "COL2A1" & all(tmp$COL2A1_check=="FAIL")){
        # if target is COL2A1 and both reps fail, sample should be listed under low_input_fail
        x <- data.frame(sample = tmp$Sample, COL2A1_Ct = tmp$ct)
        low_input_fail <- rbind(low_input_fail, x)
      }
      
      if(any(is.na(tmp$ct)) && any(!is.na(tmp$ct))){
        # Carry one Ct value forward
        und <- sum(is.na(tmp$ct))
        ind1 <- na.omit(data$Sample == s & data$Target == t)
        ct <- ifelse(und < 2, na.omit(tmp$ct), NA)
        data[ind1,]$ct <- ct
        
        if(t == "COL2A1"){
          # reprocessing required
          x <- data.frame(sample = tmp$Sample[is.na(tmp$ct)],
                          rep = tmp$rep[is.na(tmp$ct)],
                          target = t,
                          `COL2A1_check (this replicate)` = tmp$COL2A1_check[is.na(tmp$ct)])
          reprocess_needed <- rbind(reprocess_needed,x)
        } else {
          # reprocessing recommended
          x <- data.frame(sample = tmp$Sample[is.na(tmp$ct)],
                          rep = tmp$rep[is.na(tmp$ct)],
                          target = t,
                          `COL2A1_check (this replicate)` = tmp$COL2A1_check[is.na(tmp$ct)])
          reprocess <- rbind(reprocess,x)
        }
      }
    }
  }
  
  
  # Log messages
  write(ifelse(sum(low_input_fail$sample=="H2O")==2,
               paste0("Samples with insufficient input DNA (plus 2 H2O controls): ", nrow(low_input_fail)-2),
               paste0("Not all H2O controls failed to amplify COL2A1. Reprocessing strongly recommended.")),
        file = paste0(output, experimentname, "_log.txt"), append = T)
  
  if(!is_empty(reprocess_needed)){
    write(paste0("Samples that failed to amplify COL2A1 in one replicate (reprocessing needed): ", nrow(reprocess_needed)),
               file = paste0(output, experimentname, "_log.txt"), append = T)
  }
  
  if(!is_empty(reprocess)){
    write(paste0("Samples that failed to amplify one of the targets in one replicate (reprocessing recommended): ", nrow(reprocess)),
          file = paste0(output, experimentname, "_log.txt"), append = T)
  }
  
  data1 <- data %>%
    pivot_wider(id_cols = Sample,
                names_from = Target,
                values_from = ct,
                values_fn = mean) # mean is already calculated, but using "mean" as a function to suppress warning and provide only one value
  
  # data2 <- data %>%
  #   pivot_wider(id_cols = Sample,
  #               names_from = Target,
  #               values_from = ct,
  #               values_fn = sd)

  
  # Calculate input amounts (using calibration curve COL2A1)
  for (t in targets){
    varname <- paste0("input_",t)
    data1[varname] <- with(data1, 10^((data1[t]-intercept)/slope))
  }
  
  # Normalization: Divide calculated input amount by input amount calculated for COL2A1
  for (t in targets){
    varname <- paste0("ref_",t)
    calcname <- paste0("input_",t)
    data1[varname] <- data1[calcname]/data1["input_COL2A1"]
  }
  
  # Fully methylated DNA (SSS1 or gBlock) - only gBlock implemented
  gblock <- data1 %>%
    filter(Sample %in% c("gBlock", "gBLOCK", "gBlock (+)")) 
  
  # Calculate PMRs: 
  # For each ref gene divide normalized input sample by normalized input fully methylated DNA *100
  results <- data1
  
  for (t in targets2){
    varname <- t
    calcname <- paste0("ref_",t)
    div_block <- c(gblock[calcname])
    results[varname] <- (data1[calcname]/div_block)*100
  }
  
  results <- results %>%
    select(Sample, all_of(targets2))
  results[is.na(results)] <- 0 # set all NAs to 0
  
  # Calculate concentration input DNA based on COL2A mean Ct
  conc_input <- data1 %>% 
    as.data.frame() %>%
    mutate(input_conc = (data1$COL2A1 - intercept)/slope) %>%
    select(Sample, input_conc) %>% droplevels
  
  # Save results
  if(write.results==TRUE){
    results <- as.data.frame(results)
    
    data1 <- data1 %>%
      select(Sample, COL2A1, all_of(targets2)) %>%
      as.data.frame()
    
    wb <- createWorkbook()
    addWorksheet(wb, "PMR Values")
    writeDataTable(wb = wb, sheet = 1, x = results, rowNames=FALSE)
    addWorksheet(wb, "Ct")
    writeDataTable(wb = wb, sheet = 2, x = data1, rowNames=FALSE)
    addWorksheet(wb, "Input DNA conc")
    writeDataTable(wb=wb, sheet = 3, x= conc_input, rowNames=FALSE)
    addWorksheet(wb, "low DNA input")
    writeDataTable(wb=wb, sheet = 4, x= low_input_fail, rowNames=FALSE)
    
    sheet_nr=5
    
    if(!is_empty(reprocess_needed)){
      addWorksheet(wb, "Reprocessing needed")
      writeDataTable(wb = wb, sheet = sheet_nr, x = reprocess_needed, rowNames=FALSE)
      sheet_nr=sheet_nr+1
    }
    
    if(!is_empty(reprocess)){
      addWorksheet(wb, "Reprocessing recommended")
      writeDataTable(wb = wb, sheet = sheet_nr, x = reprocess, rowNames=FALSE)
    }
    
    saveWorkbook(wb, paste0(output, experimentname, ".xlsx"))
    
    cat("\nOutput file created.")
  }
  
  return(results)
  gc()
}
