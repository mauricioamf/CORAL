library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(scales)

setwd("~/Doutorado/4.Underground_metabolism/underground_metabolism/Code")

data <- read.csv("FSmatrix_SKO_glucose_TurNuP_Delta.csv", header = TRUE, sep = "\t")

#data_long <- gather(data, key = "Metabolite", value = "Basal_FS", -SubpoolKO)
data_long <- gather(data, key = "Metabolite", value = "Basal_FS")

#data_long <- data_long[data_long$Basal_FS != 0, ]
data_long <- data_long %>%
  filter(!(Basal_FS >= -0.05 & Basal_FS <= 0.05))

data_long <- data_long %>% filter(!grepl("prot_pool", Metabolite))

data_long <- data_long %>% filter(!grepl("fad_c", Metabolite))
data_long <- data_long %>% filter(!grepl("fadh2_c", Metabolite))
data_long <- data_long %>% filter(!grepl("glu__L_c", Metabolite))
data_long <- data_long %>% filter(!grepl("glyclt_c_iRN1366u", Metabolite))
data_long <- data_long %>% filter(!grepl("h2o_c", Metabolite))
data_long <- data_long %>% filter(!grepl("h2o_e", Metabolite))
data_long <- data_long %>% filter(!grepl("h2o_p", Metabolite))

#data_wide <- spread(data_long, key = Metabolite, value = Basal_FS, fill = 0)
#data_long <- gather(data_wide, key = "Metabolite", value = "Basal_FS", -SubpoolKO)

rescale_values <- function(x) {
  ifelse(x < -20, -40,
         ifelse(x > 20, 40, x))
}
data_long$Basal_FS <- rescale_values(data_long$Basal_FS)

#boxplot <- ggplot(data_long, aes(x = Metabolite, y = Basal_FS)) +
#  geom_boxplot() +
#  labs(x = "Metabolite",
#       y = "Basal flux-sum")

#boxplot + theme(axis.text.x = element_text(size = 12),
#    panel.background = element_rect(fill = NA))

ggboxplot(data_long, 
          x = "Metabolite", 
          y = "Basal_FS", 
          x.text.angle = 90,
          width = 1) + 
  scale_x_discrete(expand = expansion(mult = c(-0.5, -0.5))) +
  theme(axis.text.x = element_text(size = 10))

ggboxplot(data_long,
          x = "Metabolite",
          y = "Basal_FS",
          x.text.angle = 90,
          width = 1) +
  facet_wrap(~ Metabolite,
             nrow = 4,
             scales = "free_x") +
  scale_x_discrete(expand = expansion(mult = c(-0.5, -0.5))) +
  theme(axis.text.x = element_text(size = 6))
