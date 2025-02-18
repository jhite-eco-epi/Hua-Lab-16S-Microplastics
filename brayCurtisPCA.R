library(vegan)
library(ggplot2)
library(qiime2R)

source('./loadBaseData.R')


qza_path <- file.path(Sys.getenv("BASE_DATA_PATH"), "bray_curtis_distance_matrix.qza")
bray <- read_qza(qza_path)
bray_dist = bray$data

base_data <- loadBaseData()
# Assuming 'bray_dist' is your Bray-Curtis dissimilarity matrix
pcoa_res <- cmdscale(bray_dist, eig = TRUE, k = 2)  # 2-dimensional solution

# Compute PCoA from your Bray-Curtis dissimilarity matrix
pcoa_res <- cmdscale(bray_dist, eig = TRUE, k = 2)  # 2-dimensional solution

# Create a data frame from the ordination results
ordination_df <- as.data.frame(pcoa_res$points)
colnames(ordination_df) <- c("PCoA1", "PCoA2")

# Since row names of ordination_df contain sample IDs with the "KA" prefix,
# add them as a column for merging.
ordination_df$SampleID_full <- rownames(ordination_df)

# --------------------------------------------------
# 2. Merge Base Data (base_data) with Ordination Data
# --------------------------------------------------

# Convert the numeric Sample.ID in base_data to the full ID format (e.g., 1 -> "KA001")
base_data$SampleID_full <- as.character(sprintf("KA%03d", as.numeric(base_data$Sample.ID)))

ordination_df$SampleID_full <- as.character(rownames(ordination_df))
# Merge the datasets by the common key
merged_data <- merge(base_data, ordination_df, by = "SampleID_full", all.x = TRUE)

# Optionally, check the first few rows of the merged data
head(merged_data)

# --------------------------------------------------
# 3. Plot PCoA Results with Experimental Factors
# --------------------------------------------------

# Create a PCoA plot using ggplot.
# Here we color points by Microplastic Length and shape them by Parasite presence.
ggplot(merged_data, aes(x = PCoA1, y = PCoA2,
                        color = MPF.length..Long.Short.,
                        shape = Parasite..P.parasite.NP.non.parasite.)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(fill = MPF.length..Long.Short.), geom = "polygon", alpha = 0.2) +
  labs(title = "PCoA (Bray-Curtis) of Microbial Communities",
       x = "PCoA Axis 1",
       y = "PCoA Axis 2") +
  theme_minimal()
