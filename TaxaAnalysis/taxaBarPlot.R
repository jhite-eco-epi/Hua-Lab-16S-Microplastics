#!/usr/bin/env Rscript

# -------------------------------
# plot_community.R
#
# This script reads in a level-2 CSV file containing taxonomic counts 
# and metadata, reshapes the data into long format, computes relative abundances, 
# and creates a stacked bar plot. The x-axis represents a composite factor of the 
# treatment variables (Parasite, MPF_length, Concentration), and the stacks show 
# the relative abundance of each taxon.
#
# Make sure you have the tidyverse package installed. You can install it by running:
# install.packages("tidyverse")
# -------------------------------
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

# Identify taxonomy columns (assumed to start with "d__")
taxa_cols <- grep("^d__", names(data), value = TRUE)

# Convert from wide to long format so that each row is one sample–taxon combination
data_long <- data %>%
  pivot_longer(
    cols = all_of(taxa_cols),
    names_to = "Taxon",
    values_to = "Count"
  )

# -------------------------------
# 2. Compute relative abundance per sample
# -------------------------------
# (Assuming "index" uniquely identifies each sample)
data_long <- data_long %>%
  group_by(index) %>%
  mutate(RelAbundance = Count / sum(Count, na.rm = TRUE)) %>%
  ungroup()

# -------------------------------
# 3. Ensure treatment columns are factors
# -------------------------------
# Adjust factor levels as needed.
data_long <- data_long %>%
  mutate(
    Parasite = factor(Parasite, levels = c("NP", "P"), labels = c("No Parasite", "Parasite")),
    MPF_length = factor(MPF_length, levels = c("0", "S", "L"), , labels = c("Control", "Short", "Long")),
    Concentration = factor(Concentration)
  )

# -------------------------------
# 4. Aggregate data by treatment
# -------------------------------
# Compute the mean relative abundance for each Taxon within each combination of
# Parasite, Concentration, and MPF_length.
group_summary <- data_long %>%
  group_by(Parasite, Concentration, MPF_length, Taxon) %>%
  summarize(Counts = mean(Count, na.rm = TRUE), .groups = "drop")


# -------------------------------
# 5. Create the stacked and dodged bar plot
# -------------------------------
# Mapping:
#   • x: Concentration (so only concentration appears on the axis)
#   • y: mean_RelAbundance (the height of the bars)
#   • fill: Taxon (each bar will be internally stacked by taxon)
#   • group: combo (ensures that within each Concentration the different MPF_length bars are dodged separately)
p <- ggplot(group_summary, aes(x = Concentration, y = Counts, fill = Taxon, group = MPF_length)) +
  geom_col(position="stack") +
  geom_bar(position="dodge", stat="identity") +
  facet_grid(. ~ Parasite) +
  scale_fill_viridis_d(name = "Taxa", option = "plasma") +
  labs(x = "Microplastic Concentration (µg/L)",
       y = "Mean Relative Abundance") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        panel.border = element_rect(color = "gray", fill = NA, size = 1),
        panel.spacing.x = unit(1, "lines"))


p <- ggplot(group_summary, aes(x = Concentration, y = Counts, fill = Taxon, group = MPF_length)) +
  geom_col(position="stack") +
  facet_grid(. ~ Parasite) +
  scale_fill_viridis_d(name = "Taxa", option = "plasma") +
  labs(x = "Microplastic Concentration (µg/L)",
       y = "Mean Relative Abundance") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        panel.border = element_rect(color = "gray", fill = NA, size = 1),
        panel.spacing.x = unit(1, "lines"))

# -------------------------------
# 6. Add labels above each bar for MPF_length
# -------------------------------
p <- p + geom_text(
  data = label_data,
  aes(x = Concentration, y = total + 0.05, label = MPF_length),
  position = position_dodge2(width = 0.9, preserve = "single"),
  size = 4,
  vjust = 0
)

# -------------------------------
# 7. Output the plot
# -------------------------------
print(p)