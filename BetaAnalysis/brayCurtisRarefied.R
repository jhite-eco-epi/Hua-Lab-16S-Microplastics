# -------------------------------------------------
# PERMANOVA Analysis Using Bray-Curtis Distance Matrix
# -------------------------------------------------

# Load required libraries
library(vegan)
library(dplyr)

readRenviron(".Renviron")

# 1. Load Treatment Metadata
# Adjust the path and file name as needed.
metadata_path <- file.path(Sys.getenv("BASE_DATA_PATH"),
                           "Metadata for microbiota analysis.xlsx - PES MPF paper metadata.csv")
my_data <- read.csv(metadata_path, header = TRUE, stringsAsFactors = FALSE)

# Create a full sample ID (assuming sample IDs need formatting as "KA001", etc.)
my_data$SampleID_full <- sprintf("KA%03d", as.numeric(my_data$Sample.ID))

# Convert treatment variables to factors
my_data$MicroplasticLength <- factor(my_data$MPF.length..Long.Short., levels = c("0", "S", "L"))
my_data$MicroplasticConcentration <- factor(my_data$Concentration..µg.L., levels = c("0", "10", "40"))
my_data$Parasite <- factor(my_data$Parasite..P.parasite.NP.non.parasite., levels = c("NP", "P"))

# Create a binary variable for microplastic presence (control vs. treatment)
my_data$MicroplasticPresence <- ifelse(my_data$Concentration..µg.L. == "0", "Absent", "Present")
my_data$MicroplasticPresence <- factor(my_data$MicroplasticPresence, levels = c("Absent", "Present"))

# 2. Load Bray-Curtis Distance Matrix
# The file is expected to be a TSV with row names corresponding to SampleID_full.
bray_file <- file.path(Sys.getenv("BASE_DATA_PATH"), "rarefied_bray_curtis_distance.tsv")
bray_matrix <- read.table(bray_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Convert the matrix to a 'dist' object (required by adonis)
bray_dist <- as.dist(bray_matrix)

# 3. Match Metadata with Distance Matrix Samples
# Identify common samples between metadata and distance matrix
common_samples <- intersect(rownames(bray_matrix), my_data$SampleID_full)
my_data_subset <- my_data[my_data$SampleID_full %in% common_samples, ]

# Reorder the metadata to match the order of the distance matrix rows, if needed
my_data_subset <- my_data_subset[match(rownames(bray_matrix), my_data_subset$SampleID_full), ]

# 4. Run PERMANOVA

# Remove rows with missing values in the key variables
my_data_subset_complete <- na.omit(my_data_subset[, c("SampleID_full", "MicroplasticPresence", "Parasite")])

common_samples_complete <- intersect(rownames(bray_matrix), my_data_subset_complete$SampleID_full)
my_data_subset_complete <- my_data_subset_complete[match(common_samples_complete, my_data_subset_complete$SampleID_full), ]
bray_matrix_complete <- as.dist(bray_matrix[common_samples_complete, common_samples_complete])

# Make sure your metadata includes MicroplasticPresence and Parasite as factors
my_data_subset_complete$MicroplasticPresence <- factor(my_data_subset_complete$MicroplasticPresence, levels = c("Absent", "Present"))
my_data_subset_complete$Parasite <- factor(my_data_subset_complete$Parasite, levels = c("NP", "P"))


# Run PERMANOVA with just MicroplasticPresence and Parasite (and their interaction)
set.seed(123)  # For reproducibility of permutations
adonis_result_simple <- adonis2(bray_matrix_complete ~ MicroplasticPresence * Parasite,
                                data = my_data_subset_complete,
                                permutations = 999,
                                by = "terms")

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999

# adonis2(formula = bray_matrix_complete ~ MicroplasticPresence * Parasite, data = my_data_subset_complete, permutations = 999, by = "terms")
#                               Df SumOfSqs      R2      F Pr(>F)    
# MicroplasticPresence           1   0.6856 0.03255 3.2216  0.001 ***
# Parasite                       1   0.1836 0.00872 0.8630  0.557    
# MicroplasticPresence:Parasite  1   0.1914 0.00909 0.8993  0.516    
# Residual                      94  20.0037 0.94965                  
# Total                         97  21.0642 1.00000                  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# 1. Subset the Metadata for Treated Samples
my_data_subset_complete <- na.omit(my_data_subset[, c("SampleID_full", "MicroplasticPresence", "MicroplasticLength", "MicroplasticConcentration", "Parasite")])
# treated_data <- treated_data[complete.cases(treated_data[, c("SampleID_full", "MicroplasticLength", "MicroplasticConcentration", "Parasite")]), ]
treated_data <- subset(my_data_subset_complete, MicroplasticPresence == "Present")

# Drop any unused factor levels
treated_data$MicroplasticLength <- droplevels(treated_data$MicroplasticLength)
treated_data$MicroplasticConcentration <- droplevels(treated_data$MicroplasticConcentration)
treated_data$Parasite <- droplevels(treated_data$Parasite)

# 2. Subset the Bray-Curtis Distance Matrix for Treated Samples
# Identify common samples between the treated metadata and the distance matrix
common_treated <- intersect(rownames(bray_matrix), treated_data$SampleID_full)

# Reorder the metadata to match the distance matrix order
treated_data <- treated_data[match(common_treated, treated_data$SampleID_full), ]
treated_data <- treated_data[, c("SampleID_full", "MicroplasticLength", "MicroplasticConcentration", "Parasite")]

# Subset the distance matrix to only include the treated samples
bray_matrix_treated <- as.dist(as.matrix(bray_matrix)[common_treated, common_treated])


# 3. Run PERMANOVA on the Treated Samples
# Here we test how MicroplasticLength, MicroplasticConcentration, and Parasite (and their interactions)
# affect community composition among treated samples.
set.seed(123)  # For reproducibility
adonis_treated <- adonis2(bray_matrix_treated ~ MicroplasticLength * MicroplasticConcentration * Parasite,
                          data = treated_data,
                          permutations = 999,
                          by = "terms")
print(adonis_treated)

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999

# adonis2(formula = bray_matrix_treated ~ MicroplasticLength * MicroplasticConcentration * Parasite, data = treated_data, permutations = 999, by = "terms")
#                                                       Df SumOfSqs      R2      F Pr(>F)  
# MicroplasticLength                                     1   0.2259 0.01323 1.0494  0.353  
# MicroplasticConcentration                              1   0.2026 0.01187 0.9414  0.469  
# Parasite                                               1   0.2468 0.01445 1.1467  0.290  
# MicroplasticLength:MicroplasticConcentration           1   0.4417 0.02587 2.0524  0.030 *
# MicroplasticLength:Parasite                            1   0.2619 0.01534 1.2168  0.235  
# MicroplasticConcentration:Parasite                     1   0.2381 0.01394 1.1061  0.335  
# MicroplasticLength:MicroplasticConcentration:Parasite  1   0.1758 0.01029 0.8167  0.608  
# Residual                                              71  15.2804 0.89500                
# Total                                                 78  17.0730 1.00000                
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
