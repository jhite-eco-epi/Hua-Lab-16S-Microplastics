# Load required libraries
library(qiime2R)
library(ggplot2)
library(lme4)
library(lmerTest)
library(dplyr)
# For post-hoc Dunn's test:
library(dunn.test)

readRenviron(".Renviron")

# -------------------------------
# Load and Prepare Shannon Diversity Data
# -------------------------------

# Read the exported Shannon diversity TSV file
tsv_path <- file.path(Sys.getenv("BASE_DATA_PATH"), "rarefied_shannon.tsv")
tsv_data <- read.table(tsv_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(tsv_data)

# Rename the first column to "SampleID_full" and the second column to "shannon"
names(tsv_data)[1] <- "SampleID_full"
names(tsv_data)[2] <- "shannon"

# -------------------------------
# Load Treatment Metadata and Merge
# -------------------------------

# Load the treatments CSV file (metadata)
treatment_data <- read.csv(file.path(Sys.getenv("BASE_DATA_PATH"),
                                     "Metadata for microbiota analysis.xlsx - PES MPF paper metadata.csv"),
                           header = TRUE, stringsAsFactors = FALSE)
head(treatment_data)

# Create a matching SampleID_full in the treatment data
treatment_data$SampleID_full <- sprintf("KA%03d", as.numeric(treatment_data$Sample.ID))

# Merge treatment metadata with Shannon diversity data
merged_data <- merge(treatment_data, tsv_data, by = "SampleID_full", all.x = TRUE)
head(merged_data)

# -------------------------------
# Prepare Treatment Variables
# -------------------------------

# Convert treatment variables to factors using appropriate columns
merged_data$MicroplasticLength <- factor(merged_data$MPF.length..Long.Short., levels = c("0", "S", "L"))
merged_data$MicroplasticConcentration <- factor(merged_data$Concentration..µg.L., levels = c("0", "10", "40"))
merged_data$Parasite <- factor(merged_data$Parasite..P.parasite.NP.non.parasite., levels = c("NP", "P"))



output_plot <- ggplot(merged_data, aes(x = MicroplasticLength, y = shannon, color = MicroplasticConcentration)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_jitter(position = position_dodge(width = 0.75), alpha = 0.7) +
  facet_wrap(~ Parasite) +
  labs(title = "Shannon Diversity by Length, Concentration, & Parasite",
       x = "Microplastic Length",
       y = "Shannon Diversity") +
  theme_minimal()

output_plot

library(viridis)

output_plot <- ggplot(merged_data, aes(x = MicroplasticConcentration, y = shannon, fill = MicroplasticLength)) +
  geom_boxplot(position = position_dodge(width = 0.75), outlier.shape = NA) +
  geom_jitter(position = position_dodge(width = 0.75), alpha = 0.7, size = 2) +
  facet_grid(. ~ Parasite, 
             labeller = labeller(Parasite = c("NP" = "No Parasite", "P" = "Parasite"))) +
  scale_fill_viridis_d(name = "Microplastic Length",
                       labels = c("Control", "Short", "Long"),
                       option = "plasma") +
  labs(x = "Microplastic Concentration (µg/L)",
       y = "Shannon Diversity Index",
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

print(output_plot)


# Create a new combined treatment variable based on concentration and length.
# If concentration is "0", then label as "Control". Otherwise, label as "Treated_{Concentration}_{Length}"
merged_data$CombinedTreatment <- ifelse(merged_data$MicroplasticConcentration == "0", 
                                    "Control", 
                                    paste("Treated",
                                          merged_data$MicroplasticConcentration,
                                          merged_data$MicroplasticLength,
                                          sep = "_"))
# Convert the new variable to a factor
merged_data$CombinedTreatment <- factor(merged_data$CombinedTreatment)

library(rcompanion)

scheirerRayHare(shannon ~ CombinedTreatment * Parasite, data = merged_data)


# DV:  shannon 
# Observations:  98 
# D:  1 
# MS total:  808.5 

#                            Df Sum Sq      H p.value
# CombinedTreatment           4   2142 2.6490 0.61816
# Parasite                    1    321 0.3965 0.52888
# CombinedTreatment:Parasite  4   5621 6.9526 0.13841
# Residuals                  88  70390 
