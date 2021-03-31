# install packages
# install.packages("dplyr")
library(dplyr)
# install.packages("readr")
library(readr)

# define pipeline
# pg_df: input dataframe
# removeContaminant, removeReverse: whether rows with possible contaminants or reverse sequences, respectively, are kept (if set to false) or removed (if set to true).
# type: type of intensity to extract (possible values: "MS2", "LFQ", "Raw")
# log2t: whether intensities are transformed using log2.
# one_protein_per_row: if set to true, N proteins identified in the same row are split into N different rows.
runPipeline <- function(pg_df, removeContaminant = TRUE, removeReverse = TRUE, type = "Raw", log2t = FALSE, one_protein_per_row = FALSE) {
  # create copy of input dataframe
  out_df <- tibble(pg_df)
  
  # remove unnecessary columns
  out_df <- out_df %>%
    select(`Protein IDs`, `Gene names`,starts_with("Intensity"), starts_with("LFQ intensity"), starts_with("MS/MS count"), Reverse, `Potential contaminant`, id, `Peptide IDs`, -Intensity)
  
  # perform log2 transformation, if specified
  # 0 is transformed to -Inf
  if (log2t) {
    out_df <- out_df %>%
      mutate(across(matches("Intensity") | matches("LFQ intensity") | matches("MS/MS count"), log2))
  }
  
  # remove columns of other types of intensities
  if (type == "Raw") {
    out_df <- out_df %>%
      select(-starts_with("LFQ intensity"), -starts_with("MS/MS count"))
  } else if (type == "LFQ") {
    out_df <- out_df %>%
      select(-starts_with("Intensity"), -starts_with("MS/MS count"))
  } else if (type == "MS2") {
    out_df <- out_df %>%
      select(-starts_with("Intensity"), -starts_with("LFQ intensity"))
  } else {
    warning("type should be specified as 'Raw', 'LFQ', or 'MS2'.")
  }
  
  # remove rows with contaminants and reverses, if specified
  if (removeContaminant) {
    out_df <- out_df %>%
      filter(is.na(`Potential contaminant`))
  }
  
  if (removeReverse) {
    out_df <- out_df %>%
      filter(is.na(`Potential contaminant`))
  }
  
  # output a dataframe
  return(out_df)
}

# read protein group files
pgdata <- read_tsv("pgdata.tsv")
pgdata_processed <- runPipeline(pgdata, log2t = TRUE)
View(pgdata_processed)