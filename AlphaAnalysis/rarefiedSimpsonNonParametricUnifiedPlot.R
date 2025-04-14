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
# Load and Prepare Simpson Diversity Data
# -------------------------------

# Read the exported Simpson diversity TSV file
tsv_path <- file.path(Sys.getenv("BASE_DATA_PATH"), "rarefied_simpson.tsv")
simpson_data <- read.table(tsv_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(simpson_data)

# Rename the first column to "SampleID_full" and the second column to "simpson"
names(simpson_data)[1] <- "SampleID_full"
names(simpson_data)[2] <- "simpson"

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

# Merge treatment metadata with Simpson diversity data
merged_data <- merge(treatment_data, simpson_data, by = "SampleID_full", all.x = TRUE)
head(merged_data)

# -------------------------------
# Prepare Treatment Variables
# -------------------------------

# Convert treatment variables to factors using appropriate columns
merged_data$MicroplasticLength <- factor(merged_data$MPF.length..Long.Short., levels = c("0", "S", "L"))
merged_data$MicroplasticConcentration <- factor(merged_data$Concentration..µg.L., levels = c("0", "10", "40"))
merged_data$Parasite <- factor(merged_data$Parasite..P.parasite.NP.non.parasite., levels = c("NP", "P"))



p_simpson <- ggplot(merged_data, aes(x = MicroplasticLength, y = simpson, color = MicroplasticConcentration)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_jitter(position = position_dodge(width = 0.75), alpha = 0.7) +
  facet_wrap(~ Parasite) +
  labs(title = "Simspon Diversity by Length, Concentration, & Parasite",
       x = "Microplastic Length",
       y = "Simpson's Diversity") +
  theme_minimal()

p_simpson
