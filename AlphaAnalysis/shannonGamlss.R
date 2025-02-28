# Load required libraries
library(qiime2R)
library(ggplot2)
library(lme4)
library(gamlss)
library(lmerTest)
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

#### separate out parasite presesnce as a binary variable
# 1. Create a binary variable for microplastic presence.
my_data$MicroplasticPresence <- ifelse(my_data$MicroplasticConcentration == "0", "Absent", "Present")
my_data$MicroplasticPresence <- factor(my_data$MicroplasticPresence, levels = c("Absent", "Present"))

# 2. For the overall model (controls vs. treatments), create a pared-down data frame
model_data_overall <- my_data[, c("shannon_entropy", "MicroplasticPresence", "Parasite")]

# Check that the data has the expected dimensions and levels
table(model_data_overall$MicroplasticPresence)
table(model_data_overall$Parasite)

# 3. Fit the overall gamlss BCT model (modeling only mu, with sigma, nu, and tau as intercept-only)
model_overall <- gamlss(
  shannon_entropy ~ MicroplasticPresence * Parasite,
  family = BCT,
  data = model_data_overall,
  control = gamlss.control(n.cyc = 50)
)
summary(model_overall)
#Global Deviance:     133.3354 
#            AIC:     147.3354 
#            SBC:     165.4301 

plot(model_overall)
wp(model_overall)

model_overall_sigma <- gamlss(
  shannon_entropy ~ MicroplasticPresence * Parasite,
  sigma.formula = ~ MicroplasticPresence * Parasite,
  family = BCT,
  data = model_data_overall,
  control = gamlss.control(n.cyc = 50)
)
summary(model_overall_sigma)
#Global Deviance:     132.6276 
#            AIC:     152.6276 
#            SBC:     178.4773 

plot(model_overall_sigma)
wp(model_overall_sigma)

# 4. For the treated samples, subset the data where microplastics are present.
treated_data <- subset(my_data, MicroplasticPresence == "Present")

# Create a new data frame that only contains the necessary columns for the treated samples model.
treated_data <- treated_data[, c("shannon_entropy", "MicroplasticLength", "MicroplasticConcentration", "Parasite")]

# Drop any unused factor levels.
treated_data$MicroplasticLength <- droplevels(treated_data$MicroplasticLength)
treated_data$MicroplasticConcentration <- droplevels(treated_data$MicroplasticConcentration)
treated_data$Parasite <- droplevels(treated_data$Parasite)

# Check counts to confirm that you have the expected levels for the treatment factors.
table(treated_data$MicroplasticConcentration)
table(treated_data$MicroplasticLength)

# 5. Fit the gamlss BCT model for the treated samples (modeling only mu)
model_treated <- gamlss(
  shannon_entropy ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  family = BCT,
  data = treated_data,
  control = gamlss.control(n.cyc = 50)
)
summary(model_treated)
#Global Deviance:     99.53428 
#            AIC:     121.5343 
#            SBC:     147.5982 
plot(model_treated)
wp(model_treated)


# 5. Fit the gamlss BCT model for the treated samples (modeling only mu)
model_treated_sigma <- gamlss(
  shannon_entropy ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  sigma.formula = ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  family = BCT,
  data = treated_data,
  control = gamlss.control(n.cyc = 50)
)
summary(model_treated_sigma)
#Global Deviance:     97.85716 
#            AIC:     133.8572 
#            SBC:     176.5072
plot(model_treated_sigma)
wp(model_treated_sigma)

