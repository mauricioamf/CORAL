library(ggplot2)
library(ggpubr)
library(tidyr)

setwd("/home/mferreira/Doutorado/4.Underground_metabolism/underground_metabolism/Results/")

data <- read.csv("SKOsolutionTableDelta_glucose.csv", header = TRUE, sep = ",")

data_long <- gather(data, key = "SubpoolKO", value = "Delta", -Subpool)

ggviolin(data_long, 
          x = "SubpoolKO", 
          y = "Delta", 
          x.text.angle = 90,
          add = "jitter",
          width = 1) + 
  #scale_x_discrete(expand = expansion(mult = c(-0.5, -0.5))) +
  theme(axis.text.x = element_text(size = 10))

ggdensity(data_long, 
          x = "Delta",
          add = "mean", 
          rug = TRUE,
          color = "SubpoolKO", 
          fill = "SubpoolKO")