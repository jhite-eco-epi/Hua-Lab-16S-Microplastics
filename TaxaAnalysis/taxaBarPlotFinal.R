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
library(stringr)

# Read in the data. This assumes a comma-delimited CSV file.
# If your file is tab-delimited, replace read_csv with read_tsv.

# -------------------------------
# Load and Prepare Simpson Diversity Data
# -------------------------------

# Read the exported Simpson diversity TSV file
csv_path <- file.path(Sys.getenv("BASE_DATA_PATH"), "level-2-rarefied.csv")
data <- read_csv(csv_path)

# 2. Identify taxonomic columns (assumed to start with "d__") and reshape the data into long format
taxa_cols <- grep("^d__", names(data), value = TRUE)

data_long <- data %>%
  pivot_longer(
    cols = all_of(taxa_cols),
    names_to = "Taxon",
    values_to = "Count"
  )

# (Optional) Compute relative abundances per sample for context (not used in the percentage plot)
data_long <- data_long %>%
  group_by(index) %>%
  mutate(RelAbundance = Count / sum(Count, na.rm = TRUE)) %>%
  ungroup()

# 3. Ensure treatment columns are factors.
#    Here MPF_length is recoded using levels "C", "S", "L" which are labeled as "Control", "Short", "Long"
data_long <- data_long %>%
  mutate(
    Parasite = factor(Parasite, levels = c("NP", "P"), labels = c("No Parasite", "Parasite")),
    MPF_length = factor(MPF_length, levels = c("0", "S", "L"), labels = c("Control", "Short", "Long")),
    Concentration = factor(Concentration)
  )

# 4. Aggregate counts by treatment combination.
#    For each (Parasite, Concentration, MPF_length, Taxon), sum the counts.
group_summary <- data_long %>%
  group_by(Parasite, Concentration, MPF_length, Taxon) %>%
  summarize(Counts = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  # Create a composite variable ("combo") for unique (Concentration, MPF_length) combinations.
  mutate(combo = paste(Concentration, MPF_length, sep = "_"))

# 5. Convert counts to percentages.
#    For each unique (Parasite, Concentration, MPF_length), compute total counts and percentage per Taxon.
group_summary <- group_summary %>%
  group_by(Parasite, Concentration, MPF_length) %>%
  mutate(TotalCounts = sum(Counts),
         Pct = (Counts / TotalCounts) * 100) %>%
  ungroup()

# 6. Clean up Taxon labels:
#    Remove the prefix "d__Bacteria;p__" and replace any occurrence of "__" with "Unknown".
#group_summary <- group_summary %>%
#  mutate(Taxon_clean = str_replace(Taxon, "^d__Bacteria;p__", ""),
#         Taxon_clean = if_else(Taxon_clean == "__", "Unclassified", Taxon_clean))


group_summary <- group_summary %>%
  mutate(Taxon_clean = case_when(
    Taxon == "d__Bacteria;__" ~ "Unclassified",
    str_detect(Taxon, "^d__Bacteria;p__") ~ str_replace(Taxon, "^d__Bacteria;p__", ""),
    TRUE ~ Taxon  # If it doesn't match either case, leave it as-is.
  )) %>%
  # Also, if after removal the Taxon becomes "__", replace it with "Unclassified"
  mutate(Taxon_clean = if_else(Taxon_clean == "__", "Unclassified", Taxon_clean))


# Now reorder the Taxon_clean factor:
all_taxa <- sort(setdiff(unique(group_summary$Taxon_clean), "Unclassified"))
ordered_levels <- c("Unclassified", all_taxa)
group_summary <- group_summary %>%
  mutate(Taxon_clean = factor(Taxon_clean, levels = ordered_levels))
  
# 7. Redefine the composite variable "combo" as a factor.
#    This step is optional if you want to ensure a specific order.
#    For each Concentration, MPF_length levels will appear in the order defined in the MPF_length factor.
combo_levels <- with(data_long, 
  as.vector(t(outer(levels(Concentration), levels(MPF_length), paste, sep = "_")))
)
group_summary <- group_summary %>%
  mutate(combo = factor(paste(Concentration, MPF_length, sep = "_"), levels = combo_levels))

# 8. Prepare label data for each unique combo.
#    Since percentages sum to 100 in each bar, position the MPF_length label just above 100.
label_data <- group_summary %>%
  distinct(Parasite, Concentration, MPF_length, combo) %>%
  mutate(total_pct = 100)

# 9. Create the stacked bar plot.
#    - x-axis: composite variable "combo" (later re-labeled to show only Concentration)
#    - y-axis: Percentage (Pct)
#    - fill: Taxon_clean (cleaned taxon labels)
#    - Facet by Parasite.
p <- ggplot(group_summary, aes(x = combo, y = Pct, fill = Taxon_clean)) +
  geom_col(position = "stack", width = 0.8) +
  facet_grid(. ~ Parasite) +
  scale_x_discrete(
    name = "Microplastic Concentration (µg/L)",
    labels = function(x) sapply(strsplit(x, "_"), `[`, 1)
  ) +
  scale_y_continuous(
    name = "Percentage (%)",
    limits = c(0, 110),
    breaks = c(0, 50, 100)
  ) +
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
#     Use inherit.aes = FALSE to avoid inheriting the fill/group aesthetics.
p <- p + geom_text(
  data = label_data,
  mapping = aes(x = combo, y = total_pct + 2, label = MPF_length),
  inherit.aes = FALSE,
  size = 4,
  vjust = 0
)

# 11. Display and save the plot.
print(p)
