# -------------------------------------------------
# Complete Data Import, Processing, and NMDS Plotting
# -------------------------------------------------

# Load required libraries
library(vegan)
library(dplyr)
library(ggplot2)

readRenviron(".Renviron")

# 1. Load Treatment Metadata
metadata_path <- file.path(Sys.getenv("BASE_DATA_PATH"),
                           "Metadata for microbiota analysis.xlsx - PES MPF paper metadata.csv")
my_data <- read.csv(metadata_path, header = TRUE, stringsAsFactors = FALSE)

# Create a full sample ID (e.g., "KA001", "KA002", etc.)
my_data$SampleID_full <- sprintf("KA%03d", as.numeric(my_data$Sample.ID))

# Convert treatment variables to factors
my_data$MicroplasticLength <- factor(my_data$MPF.length..Long.Short., levels = c("0", "S", "L"))
my_data$MicroplasticConcentration <- factor(my_data$Concentration..µg.L., levels = c("0", "10", "40"))
my_data$Parasite <- factor(my_data$Parasite..P.parasite.NP.non.parasite., levels = c("NP", "P"))

# Create a binary variable for microplastic presence (control vs. treatment)
my_data$MicroplasticPresence <- ifelse(my_data$Concentration..µg.L. == "0", "Absent", "Present")
my_data$MicroplasticPresence <- factor(my_data$MicroplasticPresence, levels = c("Absent", "Present"))

# 2. Load Bray-Curtis Distance Matrix
bray_file <- file.path(Sys.getenv("BASE_DATA_PATH"), "rarefied_bray_curtis_distance.tsv")
bray_matrix <- read.table(bray_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Convert the matrix to a 'dist' object (required by vegan functions)
bray_dist <- as.dist(bray_matrix)

# 3. Match Metadata with Distance Matrix Samples
common_samples <- intersect(rownames(bray_matrix), my_data$SampleID_full)
my_data_subset <- my_data[my_data$SampleID_full %in% common_samples, ]
# Reorder metadata to match the row order of the distance matrix
my_data_subset <- my_data_subset[match(rownames(bray_matrix), my_data_subset$SampleID_full), ]

# -------------------------------------------------
# NMDS Plot 1: All Samples (MicroplasticPresence vs. Parasite)
# -------------------------------------------------

# Remove rows with missing values in key columns for complete data analysis
my_data_subset_complete <- na.omit(my_data_subset[, c("SampleID_full", "MicroplasticPresence", "Parasite")])
common_samples_complete <- intersect(rownames(bray_matrix), my_data_subset_complete$SampleID_full)
my_data_subset_complete <- my_data_subset_complete[match(common_samples_complete, my_data_subset_complete$SampleID_full), ]
bray_matrix_complete <- as.dist(bray_matrix[common_samples_complete, common_samples_complete])

# Perform NMDS on the complete Bray–Curtis distance matrix
nmds_all <- metaMDS(bray_matrix_complete, k = 2, trymax = 100)

# Extract NMDS scores and merge with metadata
scores_all <- as.data.frame(scores(nmds_all))
scores_all$SampleID <- rownames(scores_all)
scores_all <- merge(scores_all, my_data_subset_complete, by.x = "SampleID", by.y = "SampleID_full")

# Create NMDS plot: color by MicroplasticPresence and shape by Parasite
p1 <- ggplot(scores_all, aes(x = NMDS1, y = NMDS2, color = MicroplasticPresence, shape = Parasite)) +
  geom_point(size = 3) +
  labs(title = "NMDS: MicroplasticPresence vs. Parasite",
       color = "MP Presence", shape = "Parasite") +
  theme_minimal()

# Display Plot 1
print(p1)

# -------------------------------------------------
# NMDS Plot 2: Treated Samples Only (Full Model: MP Length x MP Conc x Parasite)
# -------------------------------------------------

# Remove rows with missing values for the full treatment model
my_data_subset_complete2 <- na.omit(my_data_subset[, c("SampleID_full", "MicroplasticPresence",
                                                        "MicroplasticLength", "MicroplasticConcentration", "Parasite")])

# Subset to only treated samples (MicroplasticPresence == "Present")
treated_data <- subset(my_data_subset_complete2, MicroplasticPresence == "Present")
# Drop unused factor levels
treated_data$MicroplasticLength <- droplevels(treated_data$MicroplasticLength)
treated_data$MicroplasticConcentration <- droplevels(treated_data$MicroplasticConcentration)
treated_data$Parasite <- droplevels(treated_data$Parasite)

# Identify common samples for treated data and subset the distance matrix
common_treated <- intersect(rownames(bray_matrix), treated_data$SampleID_full)
treated_data <- treated_data[match(common_treated, treated_data$SampleID_full), ]
bray_matrix_treated <- as.dist(as.matrix(bray_matrix)[common_treated, common_treated])

# Perform NMDS on the treated Bray–Curtis distance matrix
nmds_treated <- metaMDS(bray_matrix_treated, k = 2, trymax = 100)

# Extract NMDS scores and merge with treated metadata
scores_treated <- as.data.frame(scores(nmds_treated))
scores_treated$SampleID <- rownames(scores_treated)
scores_treated <- merge(scores_treated, treated_data, by.x = "SampleID", by.y = "SampleID_full")

scores_treated$MicroplasticPresence <- NULL
# Create NMDS plot with faceting: rows = MicroplasticLength, columns = MicroplasticConcentration, shape = Parasite
p2 <- ggplot(scores_treated, aes(x = NMDS1, y = NMDS2, shape = Parasite)) +
  geom_point(size = 3) +
  facet_grid(MicroplasticLength ~ MicroplasticConcentration) +
  labs(title = "NMDS: Treated Samples (MP Length x MP Conc x Parasite)",
       shape = "Parasite") +
  theme_minimal()

# Display Plot 2
print(p2)


# ----------------------------
# Plot B: Targeted Interaction Plot of MP Length x MP Concentration
# (Using mean NMDS1 score as the response)
# ----------------------------
# Summarize NMDS1 by MicroplasticLength and MicroplasticConcentration
summary_scores <- scores_treated %>%
  group_by(MicroplasticConcentration, MicroplasticLength) %>%
  summarize(mean_NMDS1 = mean(NMDS1),
            se_NMDS1 = sd(NMDS1) / sqrt(n()),
            .groups = "drop")

# Create interaction plot: x-axis = concentration, separate lines for length
p_interaction <- ggplot(summary_scores, 
                        aes(x = MicroplasticConcentration, y = mean_NMDS1,
                            group = MicroplasticLength, color = MicroplasticLength)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_NMDS1 - se_NMDS1, ymax = mean_NMDS1 + se_NMDS1), width = 0.1) +
  labs(title = "Interaction: MP Length by Concentration (Mean NMDS1)",
       x = "Microplastic Concentration (µg/L)",
       y = "Mean NMDS1",
       color = "MP Length") +
  theme_minimal()

# Display the targeted interaction plot
print(p_interaction)
