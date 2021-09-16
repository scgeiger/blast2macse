pacman::p_load(tidyverse, janitor, RColorBrewer, gridExtra, plotly)

SUMMARY.FILE <- "gene-length.tsv"
SUMMARY.PREFIX <- "gene-length"
print(SUMMARY.PREFIX)

d <- read.table(SUMMARY.FILE, sep = "\t", as.is = TRUE)
acrR_df <- subset(d, V1 == "acrR")

n <-  as_tibble(d)

plot <- ggplot(n, aes(x = V2) +
#    geom_text(aes(label=V1), hjust = 0, vjust = 0) +
    theme_light() +
    geom_boxplot() +
#    geom_point(size = 3, aes(x = (V1=="acrR"), colour = "red", alpha = 0.3))+
    xlab("Gene") +
    ylab("Length")



 n %>%
  ggplot(aes(x=V1, y=V2)) +
  geom_point(alpha=0.3) +
  geom_point(data=acrR_df, aes(x=V1, y=V2, colour = 'red', size = 3)



#plotting dataset size boxplot
#need log axes
library(ggplot2, tidyverse, plyr)
d <- read.table("ST-seq-list.txt", sep = "\t", as.is = TRUE)

freq <- count(d, "V2")

p <- ggplot(freq, aes(x = freq)) +
    geom_boxplot() +
    xlab("Distribution of n Seqs in ST") +
    scale_x_log10() +
    theme_light() 
    

