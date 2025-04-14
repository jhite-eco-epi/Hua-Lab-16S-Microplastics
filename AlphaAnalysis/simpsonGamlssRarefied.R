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

tsv_path <- file.path(Sys.getenv("BASE_DATA_PATH"), "rarefied_simpson.tsv")
# Read the Simpson data (assumes tab-delimited with a header)
simpson_data <- read.table(tsv_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(simpson_data)

# Rename the first column to "SampleID_full" and the second column to "simpson"
names(simpson_data)[1] <- "SampleID_full"
names(simpson_data)[2] <- "simpson"
# Expected output:
#   SampleID_full   simpson
#   KA001           0.8831962
#   KA011           0.9256630
#   ...

# Rename the Simpson column if needed (here it's already named "simpson")

# Merge Simpson data into your existing metadata (my_data)

my_data <- read.csv(file.path(Sys.getenv("BASE_DATA_PATH"),
                                     "Metadata for microbiota analysis.xlsx - PES MPF paper metadata.csv"),
                           header = TRUE, stringsAsFactors = FALSE)
head(my_data)

# Create a matching SampleID_full in the treatment data
my_data$SampleID_full <- sprintf("KA%03d", as.numeric(my_data$Sample.ID))

# Merge treatment metadata with Simpson diversity data
my_data <- merge(my_data, simpson_data, by = "SampleID_full", all.x = TRUE)
head(my_data)


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

# Family:  c("BCT", "Box-Cox t") 

# Call:  gamlss(formula = simpson ~ MicroplasticPresence * Parasite, family = BCT,      data = model_data_overall, control = gamlss.control(n.cyc = 50)) 

# Fitting method: RS() 

# ------------------------------------------------------------------
# Mu link function:  identity
# Mu Coefficients:
#                                        Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                            0.921008   0.008026 114.756   <2e-16 ***
# MicroplasticPresencePresent            0.001650   0.008886   0.186    0.853    
# ParasiteP                              0.010199   0.011121   0.917    0.361    
# MicroplasticPresencePresent:ParasiteP -0.018039   0.012365  -1.459    0.148    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# ------------------------------------------------------------------
# Sigma link function:  log
# Sigma Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -3.20676    0.07143   -44.9   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# ------------------------------------------------------------------
# Nu link function:  identity 
# Nu Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   20.663      1.886   10.96   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# ------------------------------------------------------------------
# Tau link function:  log 
# Tau Coefficients:
#              Estimate Std. Error  t value Pr(>|t|)    
# (Intercept) 4.307e+01  1.010e-06 42641244   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# ------------------------------------------------------------------
# No. of observations in the fit:  98 
# Degrees of Freedom for the fit:  7
#       Residual Deg. of Freedom:  91 
#                       at cycle:  16 
 
# Global Deviance:     -379.8437 
#             AIC:     -365.8437 
#             SBC:     -347.749 
# ******************************************************************



model_overall_sigma <- gamlss(
  simpson ~ MicroplasticPresence * Parasite,
  sigma.formula = ~ MicroplasticPresence * Parasite,
  family = BCT,
  data = model_data_overall,
  control = gamlss.control(n.cyc = 100)
)
summary(model_overall_sigma)


# ******************************************************************
# Family:  c("BCT", "Box-Cox t") 

# Call:  gamlss(formula = simpson ~ MicroplasticPresence * Parasite, sigma.formula = ~MicroplasticPresence *  
#     Parasite, family = BCT, data = model_data_overall, control = gamlss.control(n.cyc = 50)) 

# Fitting method: RS() 

# ------------------------------------------------------------------
# Mu link function:  identity
# Mu Coefficients:
#                                        Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                            0.890641   0.009734  91.500  < 2e-16 ***
# MicroplasticPresencePresent            0.033957   0.010413   3.261  0.00155 ** 
# ParasiteP                              0.051769   0.011187   4.627 1.18e-05 ***
# MicroplasticPresencePresent:ParasiteP -0.062901   0.012420  -5.065 2.03e-06 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# ------------------------------------------------------------------
# Sigma link function:  log
# Sigma Coefficients:
#                                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                            -2.2288     0.2357  -9.456 2.64e-15 ***
# MicroplasticPresencePresent            -1.0596     0.2609  -4.062 0.000101 ***
# ParasiteP                              -1.5906     0.3249  -4.896 4.05e-06 ***
# MicroplasticPresencePresent:ParasiteP   1.7280     0.3618   4.777 6.55e-06 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# ------------------------------------------------------------------
# Nu link function:  identity 
# Nu Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   20.537      1.555   13.21   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# ------------------------------------------------------------------
# Tau link function:  log 
# Tau Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   14.186      2.205   6.435  4.7e-09 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# ------------------------------------------------------------------
# No. of observations in the fit:  98 
# Degrees of Freedom for the fit:  10
#       Residual Deg. of Freedom:  88 
#                       at cycle:  50 
 
# Global Deviance:     -382.5917 
#             AIC:     -362.5917 
#             SBC:     -336.742 
# ******************************************************************


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
  control = gamlss.control(n.cyc = 100)
)
summary(model_treated_simpson)
plot(model_treated_simpson)
wp(model_treated_simpson)

# ******************************************************************
# Family:  c("BCT", "Box-Cox t") 

# Call:  gamlss(formula = simpson ~ MicroplasticLength * MicroplasticConcentration *  
#     Parasite, family = BCT, data = treated_data, control = gamlss.control(n.cyc = 50)) 

# Fitting method: RS() 

# ------------------------------------------------------------------
# Mu link function:  identity
# Mu Coefficients:
#                                                            Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                0.926940   0.009048 102.447   <2e-16 ***
# MicroplasticLengthL                                       -0.000281   0.010294  -0.027   0.9783    
# MicroplasticConcentration40                               -0.001357   0.010327  -0.131   0.8958    
# ParasiteP                                                 -0.018036   0.010396  -1.735   0.0873 .  
# MicroplasticLengthL:MicroplasticConcentration40           -0.010727   0.014734  -0.728   0.4691    
# MicroplasticLengthL:ParasiteP                              0.006656   0.014592   0.456   0.6497    
# MicroplasticConcentration40:ParasiteP                      0.012880   0.014686   0.877   0.3836    
# MicroplasticLengthL:MicroplasticConcentration40:ParasiteP  0.002619   0.020897   0.125   0.9006    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# ------------------------------------------------------------------
# Sigma link function:  log
# Sigma Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -3.2755     0.2401  -13.64   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# ------------------------------------------------------------------
# Nu link function:  identity 
# Nu Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   19.676      3.543   5.553 5.04e-07 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# ------------------------------------------------------------------
# Tau link function:  log 
# Tau Coefficients:
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)    12.22     228.66   0.053    0.958

# ------------------------------------------------------------------
# No. of observations in the fit:  79 
# Degrees of Freedom for the fit:  11
#       Residual Deg. of Freedom:  68 
#                       at cycle:  18 
 
# Global Deviance:     -305.0496 
#             AIC:     -283.0496 
#             SBC:     -256.9857 
# ******************************************************************




# Fit the gamlss BCT model for treated samples (modeling both μ and σ)
model_treated_simpson_sigma <- gamlss(
  simpson ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  sigma.formula = ~ MicroplasticLength * MicroplasticConcentration * Parasite,
  family = BCT,
  data = treated_data,
  control = gamlss.control(n.cyc = 100)
)
summary(model_treated_simpson_sigma)
plot(model_treated_simpson_sigma)
wp(model_treated_simpson_sigma)

# ******************************************************************
# Family:  c("BCT", "Box-Cox t") 

# Call:  gamlss(formula = simpson ~ MicroplasticLength * MicroplasticConcentration *      Parasite, sigma.formula = ~MicroplasticLength * MicroplasticConcentration *  
#     Parasite, family = BCT, data = treated_data, control = gamlss.control(n.cyc = 100)) 

# Fitting method: RS() 

# ------------------------------------------------------------------
# Mu link function:  identity
# Mu Coefficients:
#                                                            Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                0.936190   0.005603 167.079  < 2e-16 ***
# MicroplasticLengthL                                       -0.004171   0.008663  -0.482   0.6316    
# MicroplasticConcentration40                               -0.008759   0.009167  -0.955   0.3426    
# ParasiteP                                                 -0.045667   0.010522  -4.340 4.63e-05 ***
# MicroplasticLengthL:MicroplasticConcentration40           -0.043378   0.014662  -2.959   0.0042 ** 
# MicroplasticLengthL:ParasiteP                              0.035767   0.013876   2.578   0.0120 *  
# MicroplasticConcentration40:ParasiteP                      0.036909   0.014939   2.471   0.0159 *  
# MicroplasticLengthL:MicroplasticConcentration40:ParasiteP  0.010208   0.021403   0.477   0.6349    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# ------------------------------------------------------------------
# Sigma link function:  log
# Sigma Coefficients:
#                                                           Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                                -3.7953     0.2236 -16.973  < 2e-16 ***
# MicroplasticLengthL                                         0.2594     0.3162   0.820 0.414750    
# MicroplasticConcentration40                                 0.4397     0.3162   1.391 0.168720    
# ParasiteP                                                   1.1797     0.3162   3.731 0.000382 ***
# MicroplasticLengthL:MicroplasticConcentration40             0.9425     0.4472   2.108 0.038605 *  
# MicroplasticLengthL:ParasiteP                              -1.2709     0.4472  -2.842 0.005850 ** 
# MicroplasticConcentration40:ParasiteP                      -1.0211     0.4472  -2.283 0.025408 *  
# MicroplasticLengthL:MicroplasticConcentration40:ParasiteP   0.0607     0.6368   0.095 0.924333    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# ------------------------------------------------------------------
# Nu link function:  identity 
# Nu Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   20.157      1.504    13.4   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# ------------------------------------------------------------------
# Tau link function:  log 
# Tau Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 21.594786   0.001488   14510   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# ------------------------------------------------------------------
# No. of observations in the fit:  79 
# Degrees of Freedom for the fit:  18
#       Residual Deg. of Freedom:  61 
#                       at cycle:  62 
 
# Global Deviance:     -309.1864 
#             AIC:     -273.1864 
#             SBC:     -230.5363 
# ******************************************************************


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


#total interactions plot:

p_treated_simpson <- ggplot(treated_data, 
                            aes(x = MicroplasticConcentration, y = simpson, fill = MicroplasticLength)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8), alpha = 0.7) +
  geom_jitter(aes(color = MicroplasticLength), 
              position = position_dodge(width = 0.8), 
              size = 2, alpha = 0.8) +
  facet_wrap(~ Parasite) +
  labs(title = "Simpson Diversity",
       x = "Microplastic Concentration",
       y = "Simpson Diversity") +
  theme_minimal() +
  theme(legend.title = element_blank())

p_treated_simpson


