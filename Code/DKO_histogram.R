# Load necessary libraries
library(ggplot2)
library(ggpubr)
library(readxl)
library(dplyr)
library(scales)
library(readr)

#data <- read_excel("~/Doutorado/4.Underground_metabolism/underground_metabolism/Results/DKO_ratio_E1_glucose.xlsx", 
#                                   sheet = "Sheet8")

data <- read_excel("~/Doutorado/4.Underground_metabolism/underground_metabolism/Results/DKO_ratio_glucose_E1_TurNuP.xlsx", 
                   sheet = "Unaffected")

#data <- read_delim("Doutorado/4.Underground_metabolism/underground_metabolism/Results/DKO_GEM.csv", 
#                      delim = "\t", 
#                      escape_double = FALSE,
#                      col_types = cols(Ratio = col_double()),
#                      trim_ws = TRUE)
#data <- data[-c(1, 2), ]

gghistogram(data, 
            x = "Ratio", 
            fill = "HasUnd3",
            color = "HasUnd3",
            position = "dodge",
            bins = 30,
            ) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10),
                     breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000),
                     labels = c(1, 10, 100, 1000, 10000, 100000, 1000000)) +
  labs(x = "Growth rate ratio",
       y = "Log Frequency",
       fill = "Has underground reaction")
