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

tsv_path <- file.path(Sys.getenv("BASE_DATA_PATH"), "rarefied_pielou_e.tsv")
# Read the Simpson data (assumes tab-delimited with a header)
pielou_e_data <- read.table(tsv_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
head(pielou_e_data)

# Rename the first column to "SampleID_full" and the second column to "simpson"
names(pielou_e_data)[1] <- "SampleID_full"
names(pielou_e_data)[2] <- "pielou_evenness"
# Expected output:
#   SampleID_full   chao1
#   KA001           0.8831962
#   KA011           0.9256630
#   ...

# Rename the Chao1 column if needed (here it's already named "chao1")

# Merge Chao1 data into your existing metadata (my_data)

my_data <- read.csv(file.path(Sys.getenv("BASE_DATA_PATH"),
                                     "Metadata for microbiota analysis.xlsx - PES MPF paper metadata.csv"),
                           header = TRUE, stringsAsFactors = FALSE)
head(my_data)

# Create a matching SampleID_full in the treatment data
my_data$SampleID_full <- sprintf("KA%03d", as.numeric(my_data$Sample.ID))

# Merge treatment metadata with Simpson diversity data
my_data <- merge(my_data, pielou_e_data, by = "SampleID_full", all.x = TRUE)
head(my_data)


# Filter out any rows missing Pielou data
my_data <- my_data[!is.na(my_data$pielou_evenness), ]

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
#plot_data <- my_data
# Convert MicroplasticConcentration (currently a factor) to a numeric variable
# (Assuming the levels are "0", "10", "40"; these will become numeric 0, 10, 40)
plot_data$MicroplasticConcentrationNumeric <- as.numeric(as.character(plot_data$MicroplasticConcentration))

plot_data <- plot_data[, c("pielou_evenness", "MicroplasticLength", "MicroplasticConcentrationNumeric", "Parasite")]
# Create the plot
library(ggplot2)

p <- ggplot(plot_data, aes(x = MicroplasticConcentrationNumeric, y = pielou_evenness, 
                             group = Parasite, color = Parasite)) +
  geom_point(size = 3) +
  geom_line(size = 1) +
  facet_wrap(~ MicroplasticLength, nrow = 1, 
             labeller = labeller(MicroplasticLength = 
                                   c("S" = "Short Microplastics", "L" = "Long Microplastics"))) +
  scale_x_continuous(name = "Microplastic Concentration (µg/L)", breaks = c(0, 10, 40)) +
  ylab("Diversity Measure (PieLou Evenness)") +
  theme_minimal() +
  theme(strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_blank())

print(p)
