library(qiime2R)
library(ggplot2)

qza_path <- file.path(Sys.getenv("BASE_DATA_PATH"), "jaccard_distance_matrix.qza")

jaccard_distance_matrix <- read_qza(qza_path)

jaccard_artifact <- read_qza(qza_path)

# Look at the structure of the data element (the distance matrix)
str(jaccard_artifact$data)

# If it's a matrix, you might want to view the first few rows
head(jaccard_artifact$data)

qza_path <- file.path(Sys.getenv("BASE_DATA_PATH"), "shannon_vector.qza")
shannon_vector <- read_qza(qza_path)

csv_path <- file.path(Sys.getenv("BASE_DATA_PATH"), "/Metadata for microbiota analysis.xlsx - PES MPF paper metadata.csv")
my_data <- read.csv(csv_path, header = TRUE, stringsAsFactors = FALSE)

# View the first few rows of the data frame
head(my_data)


# add these data frames together!

my_data$SampleID_full <- sprintf("KA%03d", as.numeric(my_data$Sample.ID))

# Convert the shannon_vector$data (a data frame) row names to a column for merging
shannon_df <- shannon_vector$data
shannon_df$SampleID_full <- rownames(shannon_df)

# Merge the shannon_entropy field into my_data using the common key
my_data <- merge(my_data, shannon_df[, c("SampleID_full", "shannon_entropy")],
                 by = "SampleID_full", all.x = TRUE)

# Check the first few rows
head(my_data)


# Ensure MPF.length and Parasite columns are factors with appropriate labels
my_data$MPF.length..Long.Short. <- factor(my_data$MPF.length..Long.Short.,
  levels = c("0", "S", "L"),
  labels = c("None", "Short", "Long"))

my_data$Parasite..P.parasite.NP.non.parasite. <- factor(my_data$Parasite..P.parasite.NP.non.parasite.,
  levels = c("NP", "P"),
  labels = c("Non-Parasite", "Parasite"))

my_data$Concentration <- factor(my_data$Concentration..µg.L)

# main effect microplastic length
ggplot(my_data, aes(x = MPF.length..Long.Short., y = shannon_entropy)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.7) +
  labs(title = "Shannon Entropy by Microplastic Length",
       x = "Microplastic Length",
       y = "Shannon Entropy") +
  theme_minimal()

#  main effect concentration
ggplot(my_data, aes(x = Concentration, y = shannon_entropy)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.7) +
  labs(title = "Shannon Entropy by Microplastic Concentration",
       x = "Microplastic Concentration (µg/L)",
       y = "Shannon Entropy") +
  theme_minimal()

#  main effect parasite
ggplot(my_data, aes(x = Parasite..P.parasite.NP.non.parasite., y = shannon_entropy)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.7) +
  labs(title = "Shannon Entropy by Parasite Presence",
       x = "Parasite Presence",
       y = "Shannon Entropy") +
  theme_minimal()

# interaction effect length and concentration
ggplot(my_data, aes(x = MPF.length..Long.Short., y = shannon_entropy, color = Concentration)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_jitter(position = position_dodge(width = 0.75), alpha = 0.7) +
  labs(title = "Shannon Entropy by Microplastic Length & Concentration",
       x = "Microplastic Length",
       y = "Shannon Entropy") +
  theme_minimal()

# interaction effect length and parasite
ggplot(my_data, aes(x = MPF.length..Long.Short., y = shannon_entropy,
                    color = Parasite..P.parasite.NP.non.parasite.)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_jitter(position = position_dodge(width = 0.75), alpha = 0.7) +
  labs(title = "Shannon Entropy by Microplastic Length & Parasite",
       x = "Microplastic Length",
       y = "Shannon Entropy") +
  theme_minimal()

ggplot(my_data, aes(x = Concentration, y = shannon_entropy,
                    color = Parasite..P.parasite.NP.non.parasite.)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_jitter(position = position_dodge(width = 0.75), alpha = 0.7) +
  labs(title = "Shannon Entropy by Concentration & Parasite",
       x = "Microplastic Concentration (µg/L)",
       y = "Shannon Entropy") +
  theme_minimal()

# three way interaction
ggplot(my_data, aes(x = MPF.length..Long.Short., y = shannon_entropy, color = Concentration)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_jitter(position = position_dodge(width = 0.75), alpha = 0.7) +
  facet_wrap(~ Parasite..P.parasite.NP.non.parasite.) +
  labs(title = "Shannon Entropy by Length, Concentration, & Parasite",
       x = "Microplastic Length",
       y = "Shannon Entropy") +
  theme_minimal()
