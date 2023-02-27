# Calibrate COL2A1 Gblock standards on plate
# Output Plot, intercept and slope

calibGblock <- function(data){
  
  # check whether libraries already loaded
  require(dplyr)
  require(tidyverse)
  require(ggplot2)
  require(ggpubr)
  
  # initialize list for output multiple objects can be returned from function ----
  listCalib = list()
 
  # Generate standard curve
  plot = data %>%
    filter(Target == "COL2A1" & ((Sample %in% c("STD1", "STD2", "STD3","STD4") | 
                                    (Sample %in% c("STD_1", "STD_2", "STD_3", "STD_4")) | 
                                    (Sample %in% c("Std 1", "Std 2", "Std 3", "Std 4")) | 
                                    (Sample %in% c("Std_1", "Std_2", "Std_3", "Std_4")) | 
                                    (Sample %in% c("Std. 1", "Std. 2", "Std. 3", "Std. 4"))))) %>%
    mutate(x = case_when(Sample %in% c("STD1", "STD_1", "Std 1", "Std_1", "Std. 1") ~ 5, # concentration = log10(copy numbers/5uL)
                         Sample %in% c("STD2", "STD_2", "Std 2", "Std_2", "Std. 2") ~ 4,
                         Sample %in% c("STD3", "STD_3", "Std 3", "Std_3", "Std. 3") ~ 3,
                         Sample %in% c("STD4", "STD_4", "Std 4", "Std_4", "Std. 4") ~ 2)) %>%
    mutate(ct = ifelse(Cq == "Undetermined", NA, as.numeric(Cq)))
  
  fit = lm(ct ~ x, data = plot)
  anno = paste0("R2 = ", signif(summary(fit)$r.squared, 3),
                 "\ny = ", round(fit$coefficients[2],2) , "x + ", round(fit$coefficients[1],2))
  y = max(plot$ct) - 0.05*max(plot$ct)
  intercept = unname(fit$coefficients[1])
  slope = unname(fit$coefficients[2])
  R2 = summary(fit)$r.squared
  
  # Plot calibration curve
  curveplot = plot %>%
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
    xlab("Concentration \n [ log(copy number/5uL) ]") +
    ylab("Ct value")  +
    annotate("text",
             x = 4,
             y = 28,
             label = anno,
             hjust = 0,
             size = 3,
             colour = ifelse(R2 < 0.99, "red", "black"))
  
  # add plot to list
  listCalib[[1]] = curveplot
  listCalib[[2]] = intercept
  listCalib[[3]] = slope
  
  # out
  return(listCalib)
  
}