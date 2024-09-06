# Load necessary libraries
library(dplyr)
library(readxl)
library(ggplot2)

# Assuming you have loaded your data into two data frames: table1 and table2
# Let's create some sample data frames for demonstration purposes
# Replace this with your actual data loading code
table_GEM <- read_excel("DKO_GEM_growthaffected.xlsx")

table_ecGEM <- read_excel("DKO_ratio_E1_glucose.xlsx", sheet = "Affected")

# Identify common proteins
common_proteins <- intersect(table_GEM$Pair, table_ecGEM$Pair)

# Identify proteins unique to each table
unique_proteins_table_GEM <- setdiff(table_GEM$Pair, table_ecGEM$Pair)
unique_proteins_table_ecGEM <- setdiff(table_ecGEM$Pair, table_GEM$Pair)

common_proteins_ratios_table_GEM <- table_GEM[table_GEM$Pair %in% common_proteins, ]
common_proteins_ratios_table_ecGEM <- table_ecGEM[table_ecGEM$Pair %in% common_proteins, ]

common_proteins_ratios_table_GEM <- subset(common_proteins_ratios_table_GEM, select = c("Pair", "Ratio"))
common_proteins_ratios_table_ecGEM <- subset(common_proteins_ratios_table_ecGEM, select = c("Pair", "Ratio"))

common_proteins_ratios <- merge(common_proteins_ratios_table_GEM, common_proteins_ratios_table_ecGEM, by = "Pair")
colnames(common_proteins_ratios)[2] <- "Ratio_GEM"
colnames(common_proteins_ratios)[3] <- "Ratio_ecGEM"

# Plot scatterplot
ggplot(common_proteins_ratios, aes(x = Ratio_GEM, y = Ratio_ecGEM)) +
  geom_point() +
  labs(x = "Ratio GEM", y = "Ratio ecGEM")

# Calculate Pearson correlation
pearson_corr <- cor(merged_table$Outcome_1, merged_table$Outcome_2, method = "pearson")

# Calculate Spearman correlation
spearman_corr <- cor(merged_table$Outcome_1, merged_table$Outcome_2, method = "spearman")

# Print correlations
cat("Pearson correlation:", pearson_corr, "\n")
cat("Spearman correlation:", spearman_corr, "\n")

# Create a contingency matrix
contingency_matrix <- matrix(0, nrow = 2, ncol = 2,
                             dimnames = list(c("GEM", "ecGEM"),
                                             c("Common Proteins", "Unique Proteins")))

contingency_matrix["GEM", "Common Proteins"] <- length(common_proteins)
contingency_matrix["GEM", "Unique Proteins"] <- length(unique_proteins_table_GEM)
contingency_matrix["ecGEM", "Common Proteins"] <- length(common_proteins)
contingency_matrix["ecGEM", "Unique Proteins"] <- length(unique_proteins_table_ecGEM)

# Print the results
cat("Common Proteins:", common_proteins, "\n")
cat("Proteins unique to Table GEM:", unique_proteins_table_GEM, "\n")
cat("Proteins unique to Table ecGEM:", unique_proteins_table2, "\n\n")
print("Contingency Matrix:")
print(contingency_matrix)
