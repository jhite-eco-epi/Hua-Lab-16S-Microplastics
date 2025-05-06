#!/usr/bin/env Rscript
# plot_community_stacked_dodged.R
#
# This script reads in a CSV file with taxonomic counts and sample metadata,
# aggregates the data by Parasite, Microplastic Concentration, MPF_length, and Taxon,
# and then creates a stacked bar plot where each unique combination (Concentration + MPF_length)
# is represented as a separate (dodged) bar with taxa stacked within it.
# The x‑axis tick labels are relabeled to show only the Concentration,
# and a text label above each bar shows the MPF_length.

readRenviron(".Renviron")
# Load required libraries
library(tidyverse)

# Read in the data. This assumes a comma-delimited CSV file.
# If your file is tab-delimited, replace read_csv with read_tsv.

# -------------------------------
# Load and Prepare Simpson Diversity Data
# -------------------------------

# Read the exported Simpson diversity TSV file
csv_path <- file.path(Sys.getenv("BASE_DATA_PATH"), "level-2-rarefied.csv")
data <- read_csv(csv_path)

taxa_cols <- grep("^d__", names(data), value = TRUE)

data_long <- data %>%
  pivot_longer(
    cols = all_of(taxa_cols),
    names_to = "Taxon",
    values_to = "Count"
  )

# 3. (Optional) Compute relative abundances per sample (for context; not used in percentage plot)
data_long <- data_long %>%
  group_by(index) %>%
  mutate(RelAbundance = Count / sum(Count, na.rm = TRUE)) %>%
  ungroup()

# 4. Ensure your treatment columns are factors.
#    Here we set MPF_length so that levels are "C", "S", "L" with labels "Control", "Short", "Long"
data_long <- data_long %>%
  mutate(
    Parasite = factor(Parasite, levels = c("NP", "P"), labels = c("No Parasite", "Parasite")),
    MPF_length = factor(MPF_length, levels = c("0", "S", "L"), labels = c("Control", "Short", "Long")),
    Concentration = factor(Concentration)
  )

# 5. Aggregate counts by treatment combination.
#    For each (Parasite, Concentration, MPF_length, Taxon), sum the counts.
group_summary <- data_long %>%
  group_by(Parasite, Concentration, MPF_length, Taxon) %>%
  summarize(Counts = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  # Create a composite variable (combo) by concatenating Concentration and MPF_length.
  # This will be re-leveled in the next step.
  mutate(combo = paste(Concentration, MPF_length, sep = "_"))

# 6. Convert counts to percentages.
#    For each unique (Parasite, Concentration, MPF_length) group, compute the total counts,
#    then calculate the percentage for each Taxon so that percentages across taxa sum to 100.
group_summary <- group_summary %>%
  group_by(Parasite, Concentration, MPF_length) %>%
  mutate(TotalCounts = sum(Counts),
         Pct = (Counts / TotalCounts) * 100) %>%
  ungroup()

# 7. Redefine the composite variable "combo" as a factor with levels defined by the desired order.
#    This ensures that within each Concentration, the MPF_length groups appear in the order
#    specified in MPF_length's factor (i.e., "Control", "Short", "Long").
combo_levels <- with(data_long, 
  as.vector(t(outer(levels(Concentration), levels(MPF_length), paste, sep = "_")))
)
group_summary <- group_summary %>%
  mutate(combo = factor(paste(Concentration, MPF_length, sep = "_"), levels = combo_levels))

# 8. Prepare label data for each unique combo.
#    Since percentages sum to 100 in each bar, position the MPF_length label just above 100%.
label_data <- group_summary %>%
  distinct(Parasite, Concentration, MPF_length, combo) %>%
  mutate(total_pct = 100)

# 9. Create the stacked bar plot.
#    Here the x-axis uses the composite variable "combo" so that each unique (Concentration, MPF_length)
#    combination forms its own bar (with taxa stacked by Pct). Then, scale_x_discrete() is used to
#    relabel the x-axis ticks to show only the Concentration.
p <- ggplot(group_summary, aes(x = combo, y = Pct, fill = Taxon)) +
  geom_col(position = "stack", width = 0.8) +
  facet_grid(. ~ Parasite) +
  scale_x_discrete(
    name = "Microplastic Concentration (µg/L)",
    labels = function(x) sapply(strsplit(x, "_"), `[`, 1)
  ) +
  scale_y_continuous(name = "Percentage (%)", limits = c(0, 110)) +
  scale_fill_viridis_d(name = "Taxa", option = "plasma") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14),
    panel.border = element_rect(color = "gray", fill = NA, size = 1),
    panel.spacing.x = unit(1, "lines")
  )

# 10. Add text labels above each bar to indicate the MPF_length.
#      Use inherit.aes = FALSE so that only the specified mapping (x, y, label) is used.
p <- p + geom_text(
  data = label_data,
  mapping = aes(x = combo, y = total_pct + 2, label = MPF_length),
  inherit.aes = FALSE,
  size = 4,
  vjust = 0
)