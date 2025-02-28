# Load required libraries
library(qiime2R)
library(ggplot2)
library(gamlss)
library(dplyr)

# Read in environment variables from .Renviron
readRenviron(".Renviron")

# -------------------------------------------------
# 1. Load and Merge Simpson Alpha Diversity Data
# -------------------------------------------------

# Define the path to the exported Simpson TSV file
simpson_path <- file.path(Sys.getenv("BASE_DATA_PATH"), "AlphaDiversity", "exported_simpson", "alpha-diversity.tsv")

# Read the Simpson data (assumes tab-delimited with a header)
simpson_data <- read.table(simpson_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(simpson_data)
# Expected output:
#   SampleID_full   simpson
#   KA001           0.8831962
#   KA011           0.9256630
#   ...

# Rename the Simpson column if needed (here it's already named "simpson")

# Merge Simpson data into your existing metadata (my_data)

my_data <- read.csv(csv_path, header = TRUE, stringsAsFactors = FALSE)

# Create a full sample ID to match those in qza datafile
my_data$SampleID_full <- sprintf("KA%03d", as.numeric(my_data$Sample.ID))
my_data <- merge(my_data, simpson_data, by = "SampleID_full", all.x = TRUE)

# Filter out any rows missing Simpson data
my_data <- my_data[!is.na(my_data$simpson), ]

# -------------------------------------------------
# 2. Prepare Treatment Variables
# -------------------------------------------------

# Convert treatment variables to factors.
# (Replace these column names with the actual column names from your CSV.)
my_data$MicroplasticLength <- factor(my_data$MPF.length..Long.Short., levels = c("0", "S", "L"))
my_data$MicroplasticConcentration <- factor(my_data$Concentration..µg.L., levels = c("0", "10", "40"))
my_data$Parasite <- factor(my_data$Parasite..P.parasite.NP.non.parasite., levels = c("NP", "P"))

# Create a binary variable for microplastic presence (control vs. treatment)
my_data$MicroplasticPresence <- ifelse(my_data$MicroplasticConcentration == "0", "Absent", "Present")
my_data$MicroplasticPresence <- factor(my_data$MicroplasticPresence, levels = c("Absent", "Present"))

# -------------------------------------------------
# 3. Overall Model: Controls vs. Treatments
# -------------------------------------------------

# Create a pared-down data frame with only the necessary columns
model_data_overall <- my_data[, c("simpson", "MicroplasticPresence", "Parasite")]

# Check that the data has the expected levels
table(model_data_overall$MicroplasticPresence)
table(model_data_overall$Parasite)

# Fit the overall gamlss BCT model (modeling only μ; sigma, nu, and tau estimated as intercept-only)
model_overall_simpson <- gamlss(
  simpson ~ MicroplasticPresence * Parasite,
  family = BCT,
  data = model_data_overall,
  control = gamlss.control(n.cyc = 50)
)
summary(model_overall_simpson)
plot(model_overall_simpson)
wp(model_overall_simpson)

# -------------------------------------------------
# 4. Treated Samples Model: Effects Among Treatments Only
# -------------------------------------------------

# Subset the data to include only treated samples (where microplastics are present)
treated_data <- subset(my_data, MicroplasticPresence == "Present")

# Create a new data frame that includes only the necessary columns
treated_data <- treated_data[, c("simpson", "MicroplasticLength", "MicroplasticConcentration", "Parasite")]

# Drop any unused factor levels
treated_data$MicroplasticLength <- droplevels(treated_data$MicroplasticLength)
treated_data$MicroplasticConcentration <- droplevels(treated_data$MicroplasticConcentration)
treated_data$Parasite <- droplevels(treated_data$Parasite)

# Check the counts to confirm the expected levels for the treatment factors
table(treated_data$MicroplasticConcentration)
table(treated_data$MicroplasticLength)

# Fit the gamlss BCT model for treated samples (modeling only μ)
model_treated_simpson <- gamlss(
  simpson ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  family = BCT,
  data = treated_data,
  control = gamlss.control(n.cyc = 50)
)
summary(model_treated_simpson)
plot(model_treated_simpson)
wp(model_treated_simpson)



# Fit the gamlss BCT model for treated samples (modeling only μ)
model_treated_simpson_sigma <- gamlss(
  simpson ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  sigma.formula = ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  family = BCT,
  data = treated_data,
  control = gamlss.control(n.cyc = 50)
)
summary(model_treated_simpson_sigma)
plot(model_treated_simpson_sigma)
wp(model_treated_simpson_sigma)

#PLOTS

# Plot raw Simpson data by MicroplasticPresence and Parasite
p_overall_simpson <- ggplot(model_data_overall, 
                            aes(x = MicroplasticPresence, y = simpson, fill = Parasite)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = Parasite), 
              position = position_dodge(width = 0.8), 
              size = 2, alpha = 0.8) +
  labs(title = "Simpson Diversity by Microplastic Presence & Parasite",
       x = "Microplastic Presence",
       y = "Simpson Diversity") +
  theme_minimal() +
  theme(legend.title = element_blank())

p_overall_simpson


# Plot raw Simpson data for treated samples by Concentration and Length,
# faceted by Parasite status
p_treated_simpson <- ggplot(treated_data, 
                            aes(x = MicroplasticConcentration, y = simpson, fill = MicroplasticLength)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = MicroplasticLength), 
              position = position_dodge(width = 0.8), 
              size = 2, alpha = 0.8) +
  facet_wrap(~ Parasite) +
  labs(title = "Simpson Diversity in Treated Samples",
       x = "Microplastic Concentration",
       y = "Simpson Diversity") +
  theme_minimal() +
  theme(legend.title = element_blank())

p_treated_simpson


