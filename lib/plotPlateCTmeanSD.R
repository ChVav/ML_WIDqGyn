#plot mean + SD CT across all samples on plate for a specific target

plotCTPlate = function(samples,target, meanCT,SD){
  require(dplyr)
  require(tidyverse)
  require(ggplot2)
  
  theme_set(theme_minimal() +
              theme(axis.title.x=element_text(size=12),
                    axis.title.y=element_text(size=12),
                    axis.text.x=element_text(size=5),
                    axis.text.y=element_text(size=12),
                    legend.text=element_text(size=12),
                    legend.title = element_text(size=12),
                    legend.key.size= unit(6, "mm")))
  
  ### transform data frames to long format for plotting
  Ct_means <- pivot_longer(meanCT,!Sample, names_to = "Target", values_to = "Ct_mean")
  Ct_sd <- pivot_longer(SD,!Sample, names_to = "Target", values_to = "Ct_sd")
  
  ### combine dataframes, drop NAs, select samples only
  df <- full_join(Ct_means, Ct_sd) %>%
    drop_na() %>%
    filter(Sample %in% samples)
  
  ### plot output depends on number of 
  # plot processed data all targets and save
  plot <- df %>%
    filter(Target==target) %>%
    ggplot(aes(x=Sample, y=Ct_mean)) +
    geom_point()+
    geom_errorbar(aes(ymin=Ct_mean-Ct_sd, ymax=Ct_mean+Ct_sd), width=0.5) +
    #facet_wrap(~target) +
    theme(legend.title=element_blank(), axis.text.x=element_text(angle=60,hjust=1)) + 
    xlab("") + #rename axis text label
    ylab("Ct (mean)")
  
  return(plot)
  
}



