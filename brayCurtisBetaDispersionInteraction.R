# Load required libraries
library(vegan)
library(ggplot2)
library(qiime2R)

# ---------------------------
# 1. Load and Format Base Data
# ---------------------------

# Load base experimental data
source('./loadBaseData.R')
base_data <- loadBaseData()

# Ensure sample IDs match the "KA###" format
base_data$SampleID_full <- as.character(sprintf("KA%03d", as.numeric(base_data$Sample.ID)))

# Convert treatment columns to factors:
# Fiber length: order "0" (none), "S" (short), then "L" (long)
base_data$`MPF.length..Long.Short.` <- factor(base_data$`MPF.length..Long.Short.`, levels = c("0", "S", "L"))
# Fiber concentration: convert to factor (if not already)
base_data$`Concentration..µg.L.` <- factor(base_data$`Concentration..µg.L.`)
# Parasite presence: convert to factor
base_data$`Parasite..P.parasite.NP.non.parasite.` <- factor(base_data$`Parasite..P.parasite.NP.non.parasite.`)

# ---------------------------
# 2. Load Bray-Curtis Distance Matrix
# ---------------------------

qza_path <- file.path(Sys.getenv("BASE_DATA_PATH"), "bray_curtis_distance_matrix.qza")
bray <- read_qza(qza_path)
bray_dist <- bray$data

# Convert the Bray-Curtis distance matrix to a full matrix and extract sample IDs
bray_mat <- as.matrix(bray_dist)
samples_in_dist <- rownames(bray_mat)

# ---------------------------
# 3. Construct the Combined Grouping Factor
# ---------------------------

# Build a grouping factor using the three treatment variables.
# Match the sample IDs from the Bray matrix to base_data$SampleID_full.
group_vector_full <- interaction(
  base_data$`MPF.length..Long.Short.`[match(samples_in_dist, base_data$SampleID_full)],
  base_data$`Concentration..µg.L.`[match(samples_in_dist, base_data$SampleID_full)],
  base_data$`Parasite..P.parasite.NP.non.parasite.`[match(samples_in_dist, base_data$SampleID_full)],
  sep = "."
)

# Remove samples with missing treatment data
non_na_idx <- !is.na(group_vector_full)
samples_used <- samples_in_dist[non_na_idx]
group_vector_full_used <- group_vector_full[non_na_idx]

# Subset the Bray matrix to include only samples with valid treatment data and convert back to a distance object
bray_mat_subset <- bray_mat[samples_used, samples_used]
bray_dist_subset <- as.dist(bray_mat_subset)

# ---------------------------
# 4. Compute Beta Dispersion
# ---------------------------

dispersion <- betadisper(bray_dist_subset, group = group_vector_full_used)

# ---------------------------
# 5. Prepare Data for Plotting
# ---------------------------

# Create a data frame with dispersion results.
# Assume the order of dispersion$distances corresponds to samples_used.
disp_df <- data.frame(
  SampleID_full = samples_used,
  Distance = dispersion$distances,
  Group = group_vector_full_used
)

# Split the combined grouping factor into its individual components.
group_df <- do.call(rbind, strsplit(as.character(disp_df$Group), "\\."))
colnames(group_df) <- c("FiberLength", "FiberConcentration", "ParasitePresence")

# Combine these columns with disp_df
disp_df <- cbind(disp_df, group_df)

# Convert the components to factors (ensure FiberLength order is 0, S, L)
disp_df$FiberLength <- factor(disp_df$FiberLength, levels = c("0", "S", "L"))
disp_df$FiberConcentration <- factor(disp_df$FiberConcentration)
disp_df$ParasitePresence <- factor(disp_df$ParasitePresence)

# ---------------------------
# 6. Plot Beta Dispersion with Full Interactions
# ---------------------------

# Facet the plot by Parasite Presence (rows) and Fiber Concentration (columns)
ggplot(disp_df, aes(x = FiberLength, y = Distance, fill = FiberLength)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.8) +
  facet_grid(ParasitePresence ~ FiberConcentration) +
  labs(title = "Beta Diversity Dispersion by Treatment Interaction",
       x = "Fiber Length",
       y = "Distance to Centroid") +
  theme_minimal()
