library(ggplot2)
library(ggeffects)
library(qiime2R)
library(dplyr)

readRenviron(".Renviron")
# Load the QIIME 2 artifact containing Shannon entropy
qza_path <- file.path(Sys.getenv("BASE_DATA_PATH"), "shannon_vector.qza")
shannon_vector <- read_qza(qza_path)

# Load your metadata CSV file
csv_path <- file.path(Sys.getenv("BASE_DATA_PATH"), "Metadata for microbiota analysis.xlsx - PES MPF paper metadata.csv")
my_data <- read.csv(csv_path, header = TRUE, stringsAsFactors = FALSE)

# Create a full sample ID to match those in the shannon_vector artifact
my_data$SampleID_full <- sprintf("KA%03d", as.numeric(my_data$Sample.ID))

# Convert the row names of the shannon data to a column for merging
shannon_df <- shannon_vector$data
shannon_df$SampleID_full <- rownames(shannon_df)

# Merge the shannon_entropy data into the metadata using SampleID_full
my_data <- merge(my_data, shannon_df[, c("SampleID_full", "shannon_entropy")],
                 by = "SampleID_full", all.x = TRUE)

# View the first few rows of the merged data
head(my_data)

# ------------------------------
# Prepare the data for modeling
# ------------------------------

# Convert treatment variables to factors.
# Replace these column names with the actual names in your CSV if they differ.
my_data$MicroplasticLength <- factor(my_data$MPF.length..Long.Short., levels = c("0", "S", "L"))
my_data$MicroplasticConcentration <- factor(my_data$Concentration..µg.L., levels = c("0", "10", "40"))
my_data$Parasite <- factor(my_data$Parasite..P.parasite.NP.non.parasite., levels = c("NP", "P"))

# Filter out the NAs
my_data <- my_data[!is.na(my_data$shannon_entropy),]
# (Optional) If you have a grouping variable (e.g., repeated measures, subjects, or blocks), ensure it is a factor.
# Uncomment and adjust the next line if applicable:
# my_data$Subject <- factor(my_data$Subject)



#### separate out parasite presesnce as a binary variable

# 1. Create a binary variable for microplastic presence.
my_data$MicroplasticPresence <- ifelse(my_data$MicroplasticConcentration == "0", "Absent", "Present")
my_data$MicroplasticPresence <- factor(my_data$MicroplasticPresence, levels = c("Absent", "Present"))

# 2. For the overall model (controls vs. treatments), create a pared-down data frame
model_data_overall <- my_data[, c("shannon_entropy", "MicroplasticPresence", "Parasite")]



# Boxplot and jitter of raw data
p1 <- ggplot(model_data_overall, aes(x = MicroplasticPresence, y = shannon_entropy, color = Parasite)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.8) +
  labs(title = "Shannon Entropy by Microplastic Presence & Parasite",
       x = "Microplastic Presence", y = "Shannon Entropy") +
  theme_minimal()

p1


treated_data <- subset(my_data, MicroplasticPresence == "Present")

# Create a new data frame that only contains the necessary columns for the treated samples model.
treated_data <- treated_data[, c("shannon_entropy", "MicroplasticLength", "MicroplasticConcentration", "Parasite")]

# Drop any unused factor levels.
treated_data$MicroplasticLength <- droplevels(treated_data$MicroplasticLength)
treated_data$MicroplasticConcentration <- droplevels(treated_data$MicroplasticConcentration)
treated_data$Parasite <- droplevels(treated_data$Parasite)

p3 <- ggplot(treated_data, aes(x = MicroplasticConcentration, y = shannon_entropy, fill = MicroplasticLength)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = MicroplasticLength), width = 0.15, alpha = 0.8, position = position_dodge(width = 0.8)) +
  facet_wrap(~ Parasite) +
  labs(title = "Shannon Entropy in Treated Samples",
       x = "Microplastic Concentration", y = "Shannon Entropy") +
  theme_minimal()

p3

