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
# 3. Load Bray-Curtis Distance Matrix
# -------------------------------------------
bray_file <- file.path(Sys.getenv("BASE_DATA_PATH"), "rarefied_bray_curtis_distance.tsv")
bray_matrix <- read.table(bray_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
bray_dist <- as.dist(bray_matrix)

# -------------------------------------------
# 4. Match Metadata with Distance Matrix Samples
# -------------------------------------------
# Identify common sample IDs between the Bray–Curtis matrix and metadata
common_samples <- intersect(rownames(bray_matrix), my_data_subset$SampleID_full)

# Subset and reorder my_data_subset to include only those common samples, in the order they appear in the distance matrix
my_data_subset <- my_data_subset[match(common_samples, my_data_subset$SampleID_full), ]

# Also, if needed, subset the distance matrix to include only these common samples:
bray_matrix_subset <- as.dist(as.matrix(bray_matrix)[common_samples, common_samples])

# -------------------------------------------
# 5. Perform NMDS on the Bray-Curtis Distance Matrix
# -------------------------------------------
nmds_global <- metaMDS(bray_matrix_subset, k = 2, trymax = 100)
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
  labs(title = "NMDS Ordination: Parasite * CombinedTreatment",
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
permanova_results <- adonis2(bray_matrix_subset ~ Parasite * CombinedTreatment, 
                            data = metadata_filtered, 
                            permutations = 999,
                            by="terms")

# Print the PERMANOVA results
print(permanova_results)


# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999

# adonis2(formula = bray_matrix_subset ~ Parasite * CombinedTreatment, data = metadata_filtered, permutations = 999, by = "terms")
#                            Df SumOfSqs      R2      F Pr(>F)   
# Parasite                    1   0.1816 0.00862 0.8658  0.559   
# CombinedTreatment           4   1.5555 0.07385 1.8539  0.005 **
# Parasite:CombinedTreatment  4   0.8679 0.04120 1.0344  0.407   
# Residual                   88  18.4592 0.87633                 
# Total                      97  21.0642 1.00000                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(pairwiseAdonis)

# Run pairwise PERMANOVA using the Bray-Curtis distance matrix and CombinedTreatment grouping.
# Here, bray_matrix_subset is our distance matrix (aligned with the metadata) and 
# my_data_subset$CombinedTreatment is our grouping variable.
set.seed(123)  # for reproducibility of permutation tests
pairwise_results <- pairwise.adonis(bray_matrix_subset, my_data_subset$CombinedTreatment,
                                    p.adjust.m = "bonferroni")
# Print the pairwise comparisons results
print(pairwise_results)

#                           pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1  Treated_40_L vs Treated_10_S  1 0.3137012 1.3152677 0.03432751   0.198       1.00    
# 2  Treated_40_L vs Treated_10_L  1 0.3895936 1.7678931 0.04560199   0.055       0.55    
# 3  Treated_40_L vs Treated_40_S  1 0.3228202 1.3935257 0.03629585   0.161       1.00    
# 4       Treated_40_L vs Control  1 0.5176472 2.3362468 0.06094094   0.010       0.10    
# 5  Treated_10_S vs Treated_10_L  1 0.3422518 1.7030637 0.04289502   0.089       0.89    
# 6  Treated_10_S vs Treated_40_S  1 0.2497120 1.1781599 0.03007185   0.260       1.00    
# 7       Treated_10_S vs Control  1 0.6289195 3.1195509 0.07775638   0.004       0.04   .
# 8  Treated_10_L vs Treated_40_S  1 0.1162713 0.5984347 0.01550412   0.879       1.00    
# 9       Treated_10_L vs Control  1 0.5431027 2.9601723 0.07407807   0.004       0.04   .
# 10      Treated_40_S vs Control  1 0.4678302 2.4021418 0.06096475   0.010       0.10    
