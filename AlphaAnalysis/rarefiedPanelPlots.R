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

# Filter data for only short and long microplastics (exclude "0")
plot_data <- subset(my_data, MicroplasticLength %in% c("S", "L"))

# Convert MicroplasticConcentration (currently a factor) to a numeric variable
# (Assuming the levels are "0", "10", "40"; these will become numeric 0, 10, 40)
plot_data$MicroplasticConcentrationNumeric <- as.numeric(as.character(plot_data$MicroplasticConcentration))

plot_data <- plot_data[, c("simpson", "MicroplasticLength", "MicroplasticConcentrationNumeric", "Parasite")]
# Create the plot
library(ggplot2)

p <- ggplot(plot_data, aes(x = MicroplasticConcentrationNumeric, y = simpson, 
                             group = Parasite, color = Parasite)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  facet_wrap(~ MicroplasticLength, nrow = 1, 
             labeller = labeller(MicroplasticLength = 
                                   c("S" = "Short Microplastics", "L" = "Long Microplastics"))) +
  scale_x_continuous(name = "Microplastic Concentration (µg/L)", breaks = c(0, 10, 40)) +
  ylab("Diversity Measure (Simpson)") +
  theme_minimal() +
  theme(strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_blank())

print(p)
