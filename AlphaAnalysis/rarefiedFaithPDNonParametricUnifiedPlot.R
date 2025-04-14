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
tsv_path <- file.path(Sys.getenv("BASE_DATA_PATH"), "rarefied_faith_pd.tsv")
tsv_data <- read.table(tsv_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(tsv_data)

# Rename the first column to "SampleID_full" and the second column to "simpson"
names(tsv_data)[1] <- "SampleID_full"
names(tsv_data)[2] <- "faith_pd"

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
merged_data <- merge(treatment_data, tsv_data, by = "SampleID_full", all.x = TRUE)
head(merged_data)

# -------------------------------
# Prepare Treatment Variables
# -------------------------------

# Convert treatment variables to factors using appropriate columns
merged_data$MicroplasticLength <- factor(merged_data$MPF.length..Long.Short., levels = c("0", "S", "L"))
merged_data$MicroplasticConcentration <- factor(merged_data$Concentration..µg.L., levels = c("0", "10", "40"))
merged_data$Parasite <- factor(merged_data$Parasite..P.parasite.NP.non.parasite., levels = c("NP", "P"))

output_plot <- ggplot(merged_data, aes(x = MicroplasticLength, y = faith_pd, color = MicroplasticConcentration)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_jitter(position = position_dodge(width = 0.75), alpha = 0.7) +
  facet_wrap(~ Parasite) +
  labs(title = "Faith Diversity by Length, Concentration, & Parasite",
       x = "Microplastic Length",
       y = "Faith Diversity (Phylogenetic)") +
  theme_minimal()

output_plot


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

scheirerRayHare(faith_pd ~ CombinedTreatment * Parasite, data = merged_data)

# DV:  faith_pd 
# Observations:  98 
# D:  1 
# MS total:  808.5 

#                            Df Sum Sq      H p.value
# CombinedTreatment           4   3704 4.5814 0.33300
# Parasite                    1      7 0.0086 0.92592
# CombinedTreatment:Parasite  4   2870 3.5493 0.47042
# Residuals                  88  71845               
