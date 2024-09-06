library(ggplot2)
library(ggpubr)
library(ggheatmap)
library(scales)
library(reshape2)
library(tidyr)
library(dplyr)


setwd("~/Doutorado/4.Underground_metabolism/underground_metabolism/Code")

data <- read.csv("../Results/SKO_redist_Delta_Ratio_glucose.csv", header = TRUE, sep = "\t")

data_long <- gather(data, key = "Metabolite", value = "Basal_FS", -SubpoolKO)
#data_long <- data_long[data_long$Basal_FS != 0, ]
data_long <- data_long %>%
  filter(!(Basal_FS >= -0.05 & Basal_FS <= 0.05))
data_wide <- spread(data_long, key = Metabolite, value = Basal_FS, fill = 0)
data_long <- gather(data_wide, key = "Metabolite", value = "Basal_FS", -SubpoolKO)

rescale_values <- function(x) {
  ifelse(x < -10, -15,
         ifelse(x > 10, 15, x))
}
data_long$Basal_FS <- rescale_values(data_long$Basal_FS)

ggplot(data_long, aes(x = Metabolite, y = SubpoolKO, fill = Basal_FS)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "red", 
    mid = "black", 
    high = "green", 
    midpoint = 0, 
    limits = c(min(data_long$Basal_FS), max(data_long$Basal_FS)), 
    space = "Lab", 
    na.value = "black", 
    guide = "colourbar"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.text.y = element_text(size = 8)
  ) +
  labs(fill = "Basal flux-sums", x = "Metabolites", y = "Subpool KO")
