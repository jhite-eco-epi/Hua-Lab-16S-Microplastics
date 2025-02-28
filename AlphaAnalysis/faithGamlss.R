# Load required libraries
library(qiime2R)
library(ggplot2)
library(lme4)
library(gamlss)
library(lmerTest)
library(dplyr)

# Read in environment variables
readRenviron(".Renviron")

# ------------------------------
# Load and merge Faith PD data
# ------------------------------

# Load the QIIME 2 artifact containing Faith PD data
qza_path <- file.path(Sys.getenv("BASE_DATA_PATH"), "faith_pd_vector.qza")
faith_pd_vector <- read_qza(qza_path)

# Look at the data format; it should have two columns (V1 and V2)
head(faith_pd_vector$data)

# Rename the columns for clarity: V1 becomes SampleID_full, V2 is the Faith PD measure
faith_pd_df <- faith_pd_vector$data
colnames(faith_pd_df) <- c("SampleID_full", "faith_pd")

# Merge the Faith PD data into your metadata (my_data) on SampleID_full
# (Assuming my_data already has a SampleID_full column from your previous code)
my_data <- merge(my_data, faith_pd_df, by = "SampleID_full", all.x = TRUE)
head(my_data)

# Filter out any rows missing Faith PD data
my_data <- my_data[!is.na(my_data$faith_pd),]

# ------------------------------
# Prepare treatment variables (if not already done)
# ------------------------------

# Convert treatment variables to factors using the appropriate columns from your metadata.
my_data$MicroplasticLength <- factor(my_data$MPF.length..Long.Short., levels = c("0", "S", "L"))
my_data$MicroplasticConcentration <- factor(my_data$Concentration..µg.L., levels = c("0", "10", "40"))
my_data$Parasite <- factor(my_data$Parasite..P.parasite.NP.non.parasite., levels = c("NP", "P"))

# Create a binary variable for microplastic presence:
my_data$MicroplasticPresence <- ifelse(my_data$MicroplasticConcentration == "0", "Absent", "Present")
my_data$MicroplasticPresence <- factor(my_data$MicroplasticPresence, levels = c("Absent", "Present"))

# ------------------------------
# Overall Model: Controls vs. Treatments for Faith PD
# ------------------------------

# Create a pared-down data frame with only the necessary columns for the overall model
model_data_overall_pd <- my_data[, c("faith_pd", "MicroplasticPresence", "Parasite")]

# Check the counts to be sure
table(model_data_overall_pd$MicroplasticPresence)
table(model_data_overall_pd$Parasite)

# Fit the overall gamlss BCT model (modeling only μ; sigma, nu, and tau are intercept-only)
model_overall_pd <- gamlss(
  faith_pd ~ MicroplasticPresence * Parasite,
  family = BCT,
  data = model_data_overall_pd,
  control = gamlss.control(n.cyc = 50)
)
summary(model_overall_pd)
plot(model_overall_pd)
wp(model_overall_pd, ylim.all = c(-3, 3))

# ------------------------------
# Treated Samples Model: Effects Among Treated Samples Only
# ------------------------------

# Subset the data to include only samples with microplastics present
treated_data_pd <- subset(my_data, MicroplasticPresence == "Present")

# Create a pared-down data frame for treated samples
treated_data_pd <- treated_data_pd[, c("faith_pd", "MicroplasticLength", "MicroplasticConcentration", "Parasite")]

# Drop any unused factor levels
treated_data_pd$MicroplasticLength <- droplevels(treated_data_pd$MicroplasticLength)
treated_data_pd$MicroplasticConcentration <- droplevels(treated_data_pd$MicroplasticConcentration)
treated_data_pd$Parasite <- droplevels(treated_data_pd$Parasite)

# Check the counts to ensure the expected levels are present
table(treated_data_pd$MicroplasticConcentration)
table(treated_data_pd$MicroplasticLength)

# Fit the gamlss BCT model for treated samples (modeling only μ)
model_treated_pd <- gamlss(
  faith_pd ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  family = BCT,
  data = treated_data_pd,
  control = gamlss.control(n.cyc = 50)
)
summary(model_treated_pd)
plot(model_treated_pd)
wp(model_treated_pd, ylim.all = c(-3, 3))



library(ggplot2)
library(ggeffects)  # For generating predicted effects

# Boxplot + jitter of raw Faith PD data by Microplastic Presence and Parasite
p_overall_raw <- ggplot(model_data_overall_pd, aes(x = MicroplasticPresence, y = faith_pd, color = Parasite)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  labs(title = "Faith PD by Microplastic Presence & Parasite",
       x = "Microplastic Presence", y = "Faith PD") +
  theme_minimal()

p_overall_raw


# Boxplot + jitter of raw Faith PD data by Concentration and Length (with Parasite as a facet)
p_treated_raw <- ggplot(treated_data_pd, aes(x = MicroplasticConcentration, y = faith_pd, fill = MicroplasticLength)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = MicroplasticLength), size = 2,
              position = position_dodge(width = 0.8), alpha = 0.8) +
  facet_wrap(~ Parasite) +
  labs(title = "Faith PD in Treated Samples",
       x = "Microplastic Concentration", y = "Faith PD") +
  theme_minimal()

p_treated_raw
