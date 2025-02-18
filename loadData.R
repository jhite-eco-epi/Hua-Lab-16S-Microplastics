library(qiime2R)

loadData <- function (filename) {
    qza_path <- file.path(Sys.getenv("BASE_DATA_PATH"), filename)
    artifact <- read_qza(qza_path)

    head(artifact$data)

    csv_path <- file.path(Sys.getenv("BASE_DATA_PATH"), "/Metadata for microbiota analysis.xlsx - PES MPF paper metadata.csv")
    merged_data <- read.csv(csv_path, header = TRUE, stringsAsFactors = FALSE)

    # add these data frames together!
    merged_data$SampleID_full <- sprintf("KA%03d", as.numeric(merged_data$Sample.ID))

    
    artifact_df <- artifact$data
    artifact_df$SampleID_full <- rownames(artifact_df)

    # Merge the output column from qza file into merged_data using the common key
    merged_data <- merge(merged_data, artifact_df[, c("SampleID_full", "shannon_entropy")],
                 by = "SampleID_full", all.x = TRUE)
    
    return(merged_data)
}
