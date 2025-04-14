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

# Create a binary variable for microplastic presence: when concentration is "0", it's Absent
merged_data$MicroplasticPresence <- ifelse(merged_data$MicroplasticConcentration == "0", "Absent", "Present")
merged_data$MicroplasticPresence <- factor(merged_data$MicroplasticPresence, levels = c("Absent", "Present"))

# -------------------------------
# Overall Effects: Kruskal-Wallis Tests
# -------------------------------

# Test for differences in Simpson diversity by MicroplasticPresence (Control vs. Treatment)
kw_presence <- kruskal.test(simpson ~ MicroplasticPresence, data = merged_data)
print(kw_presence)

# 	Kruskal-Wallis rank sum test

# data:  simpson by MicroplasticPresence
# Kruskal-Wallis chi-squared = 1.6629, df = 1, p-value = 0.1972

# Test for differences by Parasite status
kw_parasite <- kruskal.test(simpson ~ Parasite, data = merged_data)
print(kw_parasite)

# 	Kruskal-Wallis rank sum test

# data:  simpson by Parasite
# Kruskal-Wallis chi-squared = 0.091187, df = 1, p-value = 0.7627


# -------------------------------
# Full interaction model as a single variable:
# -------------------------------

merged_data$interactionGroup <- interaction(merged_data$MicroplasticLength, 
                                              merged_data$MicroplasticConcentration, 
                                              merged_data$Parasite, sep = "_")

# Run the Kruskal-Wallis test on the combined groups
kw_result <- kruskal.test(simpson ~ interactionGroup, data = merged_data)
print(kw_result)

# 	Kruskal-Wallis rank sum test

# data:  simpson by interactionGroup
# Kruskal-Wallis chi-squared = 10.885, df = 9, p-value = 0.2837


# -------------------------------
# Effects Among Treated Samples Only
# -------------------------------

# Subset data to include only treated samples (microplastics present)
treated_data <- subset(merged_data, MicroplasticPresence == "Present")

# Drop unused factor levels
treated_data$MicroplasticLength <- droplevels(treated_data$MicroplasticLength)
treated_data$MicroplasticConcentration <- droplevels(treated_data$MicroplasticConcentration)
treated_data$Parasite <- droplevels(treated_data$Parasite)

# Check counts
table(treated_data$MicroplasticConcentration)
table(treated_data$MicroplasticLength)

# Test for differences in Simpson diversity by MicroplasticLength
kw_length <- kruskal.test(simpson ~ MicroplasticLength, data = treated_data)
print(kw_length)

# Test for differences in Simpson diversity by MicroplasticConcentration
kw_concentration <- kruskal.test(simpson ~ MicroplasticConcentration, data = treated_data)
print(kw_concentration)

# Test for differences in Simpson diversity by Parasite (within treated samples)
kw_treated_parasite <- kruskal.test(simpson ~ Parasite, data = treated_data)
print(kw_treated_parasite)

# Dunn's post-hoc test for MicroplasticConcentration among treated samples
dunn_concentration <- dunn.test(treated_data$simpson, treated_data$MicroplasticConcentration, method = "bonferroni")
print(dunn_concentration)

# Dunn's test for MicroplasticLength:
dunn_length <- dunn.test(treated_data$simpson, treated_data$MicroplasticLength, method = "bonferroni")
print(dunn_length)


# Scheirer-Ray-Hare test
library(rcompanion)

# note there is no scheirerRayHare that handles 3 factors in R, so here is our two factor test:
scheirerRayHare(simpson ~ MicroplasticPresence * Parasite, data = merged_data)


# DV:  simpson 
# Observations:  98 
# D:  1 
# MS total:  808.5 

#                               Df Sum Sq      H p.value
# MicroplasticPresence           1   1362 1.6842 0.19437
# Parasite                       1     91 0.1125 0.73735
# MicroplasticPresence:Parasite  1   3438 4.2520 0.03920
# Residuals                     94  73551               


# interaction factor combining MicroplasticPresence and Parasite
merged_data$interactionGroup <- interaction(merged_data$MicroplasticPresence, merged_data$Parasite)

# Run Dunn's test on the response variable 'simpson' across the interaction groups
dunn_result <- dunn.test(merged_data$simpson, merged_data$interactionGroup, method = "bonferroni")

# Print the results
print(dunn_result)

# Kruskal-Wallis rank sum test

# data: x and group
# Kruskal-Wallis chi-squared = 6.0274, df = 3, p-value = 0.11


#                         Comparison of x by group                            
#                                 (Bonferroni)                                  
# Col Mean-|
# Row Mean |   Absent.N   Absent.P   Present.
# ---------+---------------------------------
# Absent.P |  -1.704353
#         |     0.2649
#         |
# Present. |  -0.587846   1.601514
#         |     1.0000     0.3278
#         |
# Present. |   0.148750   2.364461   1.208179
#         |     1.0000     0.0542     0.6809

# alpha = 0.05
# Reject Ho if p <= alpha/2

# only thing that comes close to significant here is Microplastic presence x microplastic Abscence with parasite present
