# Load required libraries
library(vegan)
library(ggplot2)
library(qiime2R)

# Load base experimental data (make sure this script returns a data frame)
source('./loadBaseData.R')
base_data <- loadBaseData()

# Load the Bray-Curtis distance matrix from your QZA file
qza_path <- file.path(Sys.getenv("BASE_DATA_PATH"), "bray_curtis_distance_matrix.qza")
bray <- read_qza(qza_path)
bray_dist <- bray$data

# Ensure base_data has a full sample ID in the same "KA###" format as in the distance matrix
base_data$SampleID_full <- as.character(sprintf("KA%03d", as.numeric(base_data$Sample.ID)))

# Convert the Bray-Curtis distance matrix to a full matrix and get its sample IDs
bray_mat <- as.matrix(bray_dist)
samples_in_dist <- rownames(bray_mat)

# Create the group vector for microplastic length by matching the samples in the distance matrix
group_vector <- base_data$MPF.length..Long.Short.[match(samples_in_dist, base_data$SampleID_full)]

# Identify samples with non-missing microplastic length (group) information
non_na_idx <- !is.na(group_vector)
samples_used <- samples_in_dist[non_na_idx]
group_vector_used <- group_vector[non_na_idx]

# Subset the Bray matrix to only include samples with valid group info and convert back to a distance object
bray_mat_subset <- bray_mat[samples_used, samples_used]
bray_dist_subset <- as.dist(bray_mat_subset)

# Compute beta dispersion using the subsetted distance matrix and the microplastic length grouping vector
dispersion <- betadisper(bray_dist_subset, group = group_vector_used)

# Create a data frame with the dispersion results.
# We assume the order of dispersion$distances corresponds to samples_used.
disp_df <- data.frame(
  SampleID_full = samples_used,
  Distance = dispersion$distances,
  MPF_length = group_vector_used
)
# Ensure MPF_length is treated as a factor
disp_df$MPF_length <- factor(disp_df$MPF_length, levels = c("0", "S", "L"))

# Plot the beta dispersion using ggplot2
ggplot(disp_df, aes(x = MPF_length, y = Distance, fill = MPF_length)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.8) +
  labs(title = "Beta Diversity Dispersion by Microplastic Length",
       x = "Microplastic Length",
       y = "Distance to Centroid") +
  theme_minimal()
