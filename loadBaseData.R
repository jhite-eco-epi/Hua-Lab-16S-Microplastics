library(qiime2R)

loadBaseData <- function () {
    csv_path <- file.path(Sys.getenv("BASE_DATA_PATH"), "/Metadata for microbiota analysis.xlsx - PES MPF paper metadata.csv")
    base_data <- read.csv(csv_path, header = TRUE, stringsAsFactors = FALSE)
    return(base_data)
}
