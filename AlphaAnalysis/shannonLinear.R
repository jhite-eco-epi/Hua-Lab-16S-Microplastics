# Load required libraries
library(qiime2R)
library(ggplot2)
library(lme4)
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
# (Optional) If you have a grouping variable (e.g., repeated measures, subjects, or blocks), ensure it is a factor.
# Uncomment and adjust the next line if applicable:
# my_data$Subject <- factor(my_data$Subject)

# ------------------------------
# Fit the linear mixed model
# ------------------------------


model <- lm(shannon_entropy ~ MicroplasticLength * MicroplasticConcentration * Parasite,
              data = my_data)

# Output the summary of the model to view main and interactive effects
summary(model)


anova(model)


# ------------------------------
# Diagnostic plots
# ------------------------------

# For a quick check of residuals and model fit
par(mfrow = c(2, 2))
plot(model)


qqnorm(resid(model))
qqline(resid(model), col = "red")

# BaseModelQQPlot.pdf shows a right skewed distribution
# below we try and transform the model output to a more normal distribution


# Calculate the maximum value (ignoring NAs)
max_val <- max(my_data$shannon_entropy, na.rm = TRUE)

# Reflect the data (adding 1 to avoid log(0))
my_data$shannon_reflected <- max_val - my_data$shannon_entropy + 1

# Apply a log transformation
my_data$log_shannon_reflected <- log(my_data$shannon_reflected)

# Fit a linear model using the transformed response
model_trans <- lm(log_shannon_reflected ~ MicroplasticLength * MicroplasticConcentration * Parasite,
                  data = my_data)

# Check the QQ plot of the residuals
qqnorm(resid(model_trans))
qqline(resid(model_trans), col = "red")






model_poisson <- glm(shannon_entropy ~ MicroplasticLength * MicroplasticConcentration * Parasite,
                     family = poisson(link = "log"),
                     data = my_data)
summary(model_poisson)



library(gamlss)

model_gamlss <- gamlss(
  shannon_entropy ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  sigma.formula = ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  family = NO,  # NO: Normal distribution
  data = na.omit(my_data)
)

summary(model_gamlss)


library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = model_gamlss)  # your fitted model
plot(simulationOutput)



model_gamlss_bcpe <- gamlss(
  shannon_entropy ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  sigma.formula = ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  # Optionally, you can also specify nu and tau formulas if you suspect these parameters vary with predictors:
  # nu.formula = ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  # tau.formula = ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  family = BCPE,
  data = my_data
)


model_gamlss_bct<- gamlss(
  shannon_entropy ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  sigma.formula = ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  # Optionally, you can also specify nu and tau formulas if you suspect these parameters vary with predictors:
  # nu.formula = ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  # tau.formula = ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  family = BCT,
  data = my_data
)

summary(model_gamlss_bct)
plot(model_gamlss_bct)






library(gamlss)

# Fit a BCT model
model_bct <- gamlss(
  shannon_entropy ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  sigma.formula = ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  # Optionally, you could also try modeling nu and tau if you suspect these parameters vary with predictors:
  # nu.formula = ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  # tau.formula = ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  family = BCT,  # Using the Box-Cox t distribution
  data = my_data
)

# Check model summary to inspect estimated parameters for mu, sigma, nu, and tau
summary(model_bct)

# Plot diagnostic plots including the worm plot
plot(model_bct)

# Generate a worm plot
wp(model_bct)



# Create a new data frame with only the necessary columns
model_data <- my_data[, c("shannon_entropy", "MicroplasticLength", "MicroplasticConcentration", "Parasite")]

# Check the dimensions and number of complete cases
dim(model_data)  # should show 98 rows if filtering already kept 98 rows
sum(complete.cases(model_data))  # should ideally be 98

# Fit the model using the pared-down data frame
library(gamlss)
model_bct <- gamlss(
  # shannon_entropy ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  shannon_entropy ~ MicroplasticLength + MicroplasticConcentration + Parasite,
  # sigma.formula = ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  # Optionally, you could also try modeling nu and tau if you suspect these parameters vary with predictors:
  # nu.formula = ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  # tau.formula = ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  family = BCT,  # Using the Box-Cox t distribution
  data = model_data
)

summary(model_bct)








#### separate out parasite presesnce as a binary variable

library(gamlss)

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
plot(model_overall)
wp(model_overall)

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
plot(model_treated)
wp(model_treated)

