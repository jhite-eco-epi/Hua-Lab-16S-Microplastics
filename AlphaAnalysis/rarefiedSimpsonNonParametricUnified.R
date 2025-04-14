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

# Check the resulting levels
table(merged_data$CombinedTreatment)


library(rcompanion)

scheirerRayHare(simpson ~ CombinedTreatment * Parasite, data = merged_data)

# DV:  simpson 
# Observations:  98 
# D:  1 
# MS total:  808.5 

#                            Df Sum Sq      H p.value
# CombinedTreatment           4   2366 2.9269 0.57013
# Parasite                    1    102 0.1257 0.72295
# CombinedTreatment:Parasite  4   6360 7.8668 0.09658
# Residuals                  88  69624               

