# Load required libraries
library(vegan)
library(dplyr)
library(ggplot2)
readRenviron(".Renviron")

# -------------------------------------------
# 1. Load Treatment Metadata and Create IDs
# -------------------------------------------
metadata_path <- file.path(Sys.getenv("BASE_DATA_PATH"),
                           "Metadata for microbiota analysis.xlsx - PES MPF paper metadata.csv")
my_data <- read.csv(metadata_path, header = TRUE, stringsAsFactors = FALSE)

# Create a full sample ID (e.g., "KA001", "KA002", etc.)
my_data$SampleID_full <- sprintf("KA%03d", as.numeric(my_data$Sample.ID))

# Convert treatment variables to factors:
my_data$MicroplasticLength <- factor(my_data$MPF.length..Long.Short., levels = c("0", "S", "L"))
my_data$MicroplasticConcentration <- factor(my_data$Concentration..µg.L., levels = c("0", "10", "40"))
my_data$Parasite <- factor(my_data$Parasite..P.parasite.NP.non.parasite., levels = c("NP", "P"))

# -------------------------------------------
# 2. Create Combined Treatment Variable
# -------------------------------------------
# For samples with 0 concentration, assign "Control". For others, combine concentration and length.
my_data$CombinedTreatment <- ifelse(my_data$MicroplasticConcentration == "0",
                                    "Control",
                                    paste("Treated",
                                          my_data$MicroplasticConcentration,
                                          my_data$MicroplasticLength,
                                          sep = "_"))
my_data$CombinedTreatment <- factor(my_data$CombinedTreatment)
table(my_data$CombinedTreatment)  # Check the levels

# -------------------------------------------
# 3. Load Jaccard Distance Matrix
# -------------------------------------------
jaccard_file <- file.path(Sys.getenv("BASE_DATA_PATH"), "rarefied_jaccard_matrix.tsv")
jaccard_matrix <- read.table(jaccard_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
jaccard_dist <- as.dist(jaccard_matrix)

# -------------------------------------------
# 4. Match Metadata with Distance Matrix Samples
# -------------------------------------------
# Identify common sample IDs between the Jaccard matrix and metadata
common_samples <- intersect(rownames(jaccard_matrix), my_data$SampleID_full)

# Subset metadata to include only common samples, and reorder to match the distance matrix
my_data_subset <- my_data[my_data$SampleID_full %in% common_samples, ]
my_data_subset <- my_data_subset[match(common_samples, my_data_subset$SampleID_full), ]

# Also, subset the distance matrix to include only these common samples:
jaccard_matrix_subset <- as.dist(as.matrix(jaccard_matrix)[common_samples, common_samples])

# -------------------------------------------
# 5. Perform NMDS on the Jaccard Distance Matrix
# -------------------------------------------
nmds_global <- metaMDS(jaccard_matrix_subset, k = 2, trymax = 100)
scores_global <- as.data.frame(scores(nmds_global))
scores_global$SampleID <- rownames(scores_global)
scores_global <- merge(scores_global, my_data_subset, by.x = "SampleID", by.y = "SampleID_full")

# -------------------------------------------
# 6. Create NMDS Plot for the Full Model
# -------------------------------------------

library(viridis)

# Create the boxplot of NMDS1 using MicroplasticConcentration as the x-axis variable
# Facet by Parasite (rows) and MicroplasticLength (columns)
beta_boxplot <- ggplot(scores_global, aes(x = MicroplasticConcentration, y = NMDS1, fill = MicroplasticLength)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_jitter(position = position_dodge(width = 0.75), alpha = 0.7, size = 2) +
  facet_grid(. ~ Parasite, 
             labeller = labeller(
               Parasite = c("NP" = "No Parasite", "P" = "Parasite"),
               MicroplasticLength = c("0" = "Control", "S" = "Short", "L" = "Long")
             )) +
  scale_fill_viridis_d(name = "Microplastic Length",
                       labels = c("Control", "Short", "Long"),
                       option = "plasma") +
  labs(x = "Microplastic Concentration (µg/L)",
       y = "Ordination Score (NMDS 1)",
       title = NULL) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # Removing default legend title since we're setting it via scale_fill_viridis_d
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        panel.border = element_rect(color = "gray", fill = NA, size = 1),
        panel.spacing.x = unit(1, "lines"))

print(beta_boxplot)



# Mapping: CombinedTreatment (color) and Parasite (shape)
p_global <- ggplot(scores_global, aes(x = NMDS1, y = NMDS2, 
                                        color = CombinedTreatment, 
                                        shape = Parasite)) +
  geom_point(size = 3) +
  stat_ellipse(aes(group = interaction(CombinedTreatment, Parasite)), 
               type = "t", size = 1, linetype = "dashed") +
  labs(title = NULL,
       color = "Combined Treatment",
       shape = "Parasite Status") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
  
print(p_global)

metadata_filtered <- my_data_subset[, c("SampleID_full", "Parasite", "CombinedTreatment")]

set.seed(123)  # Set seed for reproducibility
permanova_results <- adonis2(jaccard_matrix_subset ~ Parasite * CombinedTreatment, 
                            data = metadata_filtered, 
                            permutations = 999,
                            by="terms")

# Print the PERMANOVA results
print(permanova_results)

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999

# adonis2(formula = jaccard_matrix_subset ~ Parasite * CombinedTreatment, data = metadata_filtered, permutations = 999, by = "terms")
#                            Df SumOfSqs      R2      F Pr(>F)    
# Parasite                    1   0.1729 0.00989 0.9874  0.435    
# CombinedTreatment           4   1.0629 0.06082 1.5179  0.001 ***
# Parasite:CombinedTreatment  4   0.8332 0.04768 1.1898  0.066 .  
# Residual                   88  15.4065 0.88161                  
# Total                      97  17.4755 1.00000                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(pairwiseAdonis)

# Run pairwise PERMANOVA using the Jaccard distance matrix and CombinedTreatment grouping.
set.seed(123)  # for reproducibility of permutation tests
pairwise_results <- pairwise.adonis(jaccard_matrix_subset, my_data_subset$CombinedTreatment,
                                    p.adjust.m = "bonferroni")
# Print the pairwise comparisons results
print(pairwise_results) 

#                           pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted
# 1  Treated_40_L vs Treated_10_S  1 0.2779466 1.4450334 0.03758700   0.050       0.50
# 2  Treated_40_L vs Treated_10_L  1 0.2557469 1.3942774 0.03631472   0.042       0.42
# 3  Treated_40_L vs Treated_40_S  1 0.2744210 1.4886731 0.03867821   0.041       0.41
# 4       Treated_40_L vs Control  1 0.2870397 1.5732818 0.04187236   0.021       0.21
# 5  Treated_10_S vs Treated_10_L  1 0.2861440 1.6343506 0.04123571   0.026       0.26
# 6  Treated_10_S vs Treated_40_S  1 0.3064887 1.7417069 0.04382567   0.017       0.17
# 7       Treated_10_S vs Control  1 0.3485946 2.0045419 0.05139252   0.004       0.04
# 8  Treated_10_L vs Treated_40_S  1 0.1659226 0.9918557 0.02543751   0.448       1.00
# 9       Treated_10_L vs Control  1 0.2099132 1.2723372 0.03324430   0.092       0.92
# 10      Treated_40_S vs Control  1 0.2367879 1.4273311 0.03714364   0.046       0.46
