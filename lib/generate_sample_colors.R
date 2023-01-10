# Making pretty rainbow colours for samples :)
library(grDevices)
col <- rainbow(n = 90)
samples <- c()
for (i in 1:length(col)){
  col[i] <- paste('\"RGB(', col2rgb(col[i])[1], ",", col2rgb(col[i])[2], ",", col2rgb(col[i])[3], ')\"', sep = "")
  samples <- c(samples, paste("Sample", formatC(i, width=2, flag="0"), sep="_"))
}

df <- data.frame(Sample.name=samples, colors=col)

write.csv(df, file = "colors.csv",
          quote = TRUE,
          row.names = FALSE)
