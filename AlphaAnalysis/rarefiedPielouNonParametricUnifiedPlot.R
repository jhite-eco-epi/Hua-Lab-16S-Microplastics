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
tsv_path <- file.path(Sys.getenv("BASE_DATA_PATH"), "rarefied_pielou_e.tsv")
tsv_data <- read.table(tsv_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(tsv_data)

# Rename the first column to "SampleID_full" and the second column to "simpson"
names(tsv_data)[1] <- "SampleID_full"
names(tsv_data)[2] <- "pielou_evenness"

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



output_plot <- ggplot(merged_data, aes(x = MicroplasticLength, y = pielou_evenness, color = MicroplasticConcentration)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_jitter(position = position_dodge(width = 0.75), alpha = 0.7) +
  facet_wrap(~ Parasite) +
  labs(title = "Evenness by Length, Concentration, & Parasite",
       x = "Microplastic Length",
       y = "Pielou Diversity (Evenness)") +
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

scheirerRayHare(pielou_evenness ~ CombinedTreatment * Parasite, data = merged_data)


# DV:  pielou_evenness 
# Observations:  98 
# D:  1 
# MS total:  808.5 

#                            Df Sum Sq      H p.value
# CombinedTreatment           4   1388 1.7167 0.78769
# Parasite                    1    195 0.2413 0.62323
# CombinedTreatment:Parasite  4   6494 8.0328 0.09039
# Residuals                  88  70378            
