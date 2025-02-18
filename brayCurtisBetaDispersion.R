# Load required libraries
library(vegan)
library(ggplot2)
library(qiime2R)

# ---------------------------
# 1. Load and Format Base Data
# ---------------------------
source('./loadBaseData.R')
base_data <- loadBaseData()

# Ensure sample IDs match the "KA###" format
base_data$SampleID_full <- as.character(sprintf("KA%03d", as.numeric(base_data$Sample.ID)))

# ---------------------------
# 2. Load the Bray-Curtis Distance Matrix
# ---------------------------
qza_path <- file.path(Sys.getenv("BASE_DATA_PATH"), "bray_curtis_distance_matrix.qza")
bray <- read_qza(qza_path)
bray_dist <- bray$data

# Convert the Bray-Curtis distance matrix to a full matrix and extract sample IDs
bray_mat <- as.matrix(bray_dist)
samples_in_dist <- rownames(bray_mat)

# ============================
# Plot 1: Beta Dispersion by Microplastic Length
# ============================

# Create the group vector for microplastic length by matching the samples in the distance matrix
group_vector_length <- base_data$`MPF.length..Long.Short.`[match(samples_in_dist, base_data$SampleID_full)]

# Identify samples with non-missing microplastic length info
non_na_idx_length <- !is.na(group_vector_length)
samples_used_length <- samples_in_dist[non_na_idx_length]
group_vector_length_used <- group_vector_length[non_na_idx_length]

# Subset the Bray matrix for these samples and convert back to a distance object
bray_mat_subset_length <- bray_mat[samples_used_length, samples_used_length]
bray_dist_subset_length <- as.dist(bray_mat_subset_length)

# Compute beta dispersion for microplastic length
dispersion_length <- betadisper(bray_dist_subset_length, group = group_vector_length_used)

# Create a data frame for plotting
disp_df_length <- data.frame(
  SampleID_full = samples_used_length,
  Distance = dispersion_length$distances,
  MPF_length = factor(group_vector_length_used, levels = c("0", "S", "L"))
)

# Plot beta dispersion for microplastic length
ggplot(disp_df_length, aes(x = MPF_length, y = Distance, fill = MPF_length)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.8) +
  labs(title = "Beta Diversity Dispersion by Microplastic Length",
       x = "Microplastic Length",
       y = "Distance to Centroid") +
  theme_minimal()

# ============================
# Plot 2: Beta Dispersion by Fiber Concentration
# ============================

# Create the group vector for fiber concentration using the column "Concentration..µg.L."
group_vector_conc <- base_data$`Concentration..µg.L.`[match(samples_in_dist, base_data$SampleID_full)]

# Identify samples with non-missing concentration info
non_na_idx_conc <- !is.na(group_vector_conc)
samples_used_conc <- samples_in_dist[non_na_idx_conc]
group_vector_conc_used <- group_vector_conc[non_na_idx_conc]

# Subset the Bray matrix for these samples and convert back to a distance object
bray_mat_subset_conc <- bray_mat[samples_used_conc, samples_used_conc]
bray_dist_subset_conc <- as.dist(bray_mat_subset_conc)

# Compute beta dispersion for fiber concentration
dispersion_conc <- betadisper(bray_dist_subset_conc, group = group_vector_conc_used)

# Create a data frame for plotting
disp_df_conc <- data.frame(
  SampleID_full = samples_used_conc,
  Distance = dispersion_conc$distances,
  Concentration = factor(group_vector_conc_used)  # factor levels will be in natural order
)

# Plot beta dispersion for fiber concentration
ggplot(disp_df_conc, aes(x = Concentration, y = Distance, fill = Concentration)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.8) +
  labs(title = "Beta Diversity Dispersion by Fiber Concentration",
       x = "Fiber Concentration (µg/L)",
       y = "Distance to Centroid") +
  theme_minimal()

# ============================
# Plot 3: Beta Dispersion by Parasite Presence
# ============================

# Create the group vector for parasite presence using the column "Parasite..P.parasite.NP.non.parasite."
group_vector_parasite <- base_data$`Parasite..P.parasite.NP.non.parasite.`[match(samples_in_dist, base_data$SampleID_full)]

# Identify samples with non-missing parasite info
non_na_idx_parasite <- !is.na(group_vector_parasite)
samples_used_parasite <- samples_in_dist[non_na_idx_parasite]
group_vector_parasite_used <- group_vector_parasite[non_na_idx_parasite]

# Subset the Bray matrix for these samples and convert back to a distance object
bray_mat_subset_parasite <- bray_mat[samples_used_parasite, samples_used_parasite]
bray_dist_subset_parasite <- as.dist(bray_mat_subset_parasite)

# Compute beta dispersion for parasite presence
dispersion_parasite <- betadisper(bray_dist_subset_parasite, group = group_vector_parasite_used)

# Create a data frame for plotting
disp_df_parasite <- data.frame(
  SampleID_full = samples_used_parasite,
  Distance = dispersion_parasite$distances,
  Parasite = factor(group_vector_parasite_used)
)

# Plot beta dispersion for parasite presence
ggplot(disp_df_parasite, aes(x = Parasite, y = Distance, fill = Parasite)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.8) +
  labs(title = "Beta Diversity Dispersion by Parasite Presence",
       x = "Parasite Presence",
       y = "Distance to Centroid") +
  theme_minimal()
