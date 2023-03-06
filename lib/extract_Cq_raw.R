# Author: Charlotte Vavourakis
#

# function to load CSV and extract raw Cq values for each target

extractCqRaw <- function(input_file,output_file){
  require(dplyr)
  require(tidyverse)
  
  # read in data
  data <- read.table(file = input_file,
                     skip = 22,
                     sep = ",",
                     header = TRUE)
  
  # kick out standards, controls
  data <- data %>% 
    filter(!Sample %in% c("STD1", "STD2", "STD3","STD4",
                                         "STD_1", "STD_2", "STD_3", "STD_4",
                                         "Std 1", "Std 2", "Std 3", "Std 4",
                                         "Std_1", "Std_2", "Std_3", "Std_4",
                                         "Std. 1", "Std. 2", "Std. 3", "Std. 4",
                                         "NTC_H2O","NTC", "H2O",
                                         "gBlock")) %>%
    droplevels()
  
  # collect target and sample overviews
  targets <- unique(data$Target)
  samples <- unique(data$Sample)
  
  # make a wide df per target
  data <- data %>% 
    select(-Amp.Status) %>%
    droplevels() %>%
    pivot_wider(names_from="Target",values_from="Cq", values_fn=list)
  
  # split columns for first target
  target <- targets[1]
  df <- data[target]
  df <- as.data.frame(lapply(df, function(x) t(do.call(cbind, x))))
  rownames(df) <- data$Sample
  df <- df %>% rownames_to_column(var="Sample")
  
  # repeat for remaining targets and combine result
  if (length(targets) > 1) {
    for (i in 2:length(targets)){
      target <- targets[i]
      df2 <- data[target]
      df2 <- as.data.frame(lapply(df2, function(x) t(do.call(cbind,x))))
      rownames(df2) <- data$Sample
      df2 <- df2 %>% rownames_to_column(var="Sample")
      df <- full_join(df,df2)
    }
  }
  
  # save result
  write.csv(df,file=output_file)
  
  return(df)
  
}