# install packages
library(tidyverse)
# install.packages("dplyr")
library(dplyr)
# install.packages("readr")
library(readr)
#install.packages("fuzzyjoin")
library(fuzzyjoin)

# define pipeline
# pg_df: input dataframe
# removeContaminant, removeReverse: whether rows with possible contaminants or reverse sequences, respectively, are kept (if set to false) or removed (if set to true).
# type: type of intensity to extract (possible values: "MS2", "LFQ", "Raw")
# log2t: whether intensities are transformed using log2.
# one_protein_per_row: if set to true, N proteins identified in the same row are split into N different rows.
runMQ <- function(pg_df, removeContaminant = TRUE, removeReverse = TRUE, type = "Raw", log2t = FALSE, one_protein_per_row = FALSE) {
  # create copy of input dataframe
  out_df <- tibble(pg_df)
  
  # remove unnecessary columns
  out_df <- out_df %>%
    select(`Protein IDs`, `Gene names`, starts_with("Intensity"), starts_with("LFQ intensity"), starts_with("MS/MS count"), Reverse, `Potential contaminant`, id, `Peptide IDs`, -Intensity)
  
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
  
  # split each row into more rows depending on amount of identified proteins, if one_protein_per_row = TRUE
  if (one_protein_per_row) {
    out_df <- out_df %>%
      separate_rows(`Protein IDs`, sep = ';')
  }
  
  # output a dataframe
  return(out_df)
}

matrix_filter <- function (pg_df, mf_df, sample.id.col='Sample.ID', method='percentage', filterGroup='DIAGNOSIS',percentageNA = 0.5, maxNA=0, intensity_type="Raw", log2t) {
  # throw error if percentageNA < 0 or > 100
  if (percentageNA < 0 || percentageNA > 100) {
    stop("percentageNA must be between 0 and 100.")
  }
  
  # change -Inf to NA
  is.na(pg_df) <- sapply(pg_df, is.infinite)
  
  proteingroups <- tibble(pg_df)
  
  # convert protein groups tibble to long format
  if (intensity_type == "Raw") {
    proteingroups <- proteingroups %>%
      gather("Sample", "Intensity", starts_with("Intensity"), starts_with("LFQ intensity"), starts_with("MS/MS count"))
  } else if (intensity_type == "LFQ") {
    proteingroups <- proteingroups %>%
      gather("Sample", "LFQ intensity", starts_with("Intensity"), starts_with("LFQ intensity"), starts_with("MS/MS count"))
  } else if (intensity_type == "MS2") {
    proteingroups <- proteingroups %>%
      gather("Sample", "MS/MS count", starts_with("Intensity"), starts_with("LFQ intensity"), starts_with("MS/MS count"))
  } else {
    warning("type should be specified as 'Raw', 'LFQ', or 'MS2'.")
  }
  
  # add regex column to manifest
  manifest_local <- mf_df %>%
    mutate(Sample.ID.regex = paste("(Intensity|LFQ intensity|MS/MS count)\\s", Sample.ID, ".*", sep=""))
  
  gsub(".", "[.]", manifest_local$Sample.ID.regex)
  
  # use left join to add patient diagnosis to sample ID
  combined_df <- regex_left_join(proteingroups, manifest_local, by=c(Sample = "Sample.ID.regex"))[-ncol(proteingroups) - ncol(manifest_local)] %>%
    select(-Sample.ID) %>%
    distinct()
  
  # create NARate column names for every unique diagnosis
  diagnoses <- unique(combined_df$DIAGNOSIS)
  
  # filter values that do not meet conditions
  # first, create columns of NA counts and total diagnoses for each protein group and diagnosis
  combined_df <- combined_df %>%
    mutate(MinNARate = 1, MaxNAcount = 0)
  
  if (log2t) {
    for (i in diagnoses) {
      combined_df <- combined_df %>%
        group_by(`Protein IDs`) %>%
        mutate("NACount.{{i}}" := sum(is.na(Intensity) & DIAGNOSIS == i), "NARate.{{i}}" := sum(is.na(Intensity) & DIAGNOSIS == i)/sum(DIAGNOSIS == i), MinNARate = min(MinNARate, sum(is.na(Intensity) & DIAGNOSIS == i)/sum(DIAGNOSIS == i)), MaxNAcount = max(MaxNAcount, sum(is.na(Intensity) & DIAGNOSIS == i)))
    }
  } else {
    for (i in diagnoses) {
      combined_df <- combined_df %>%
        group_by(`Protein IDs`) %>%
        mutate("NACount.{{i}}" := sum(Intensity == 0 & DIAGNOSIS == i), "NARate.{{i}}" := sum(is.na(Intensity) & DIAGNOSIS == i)/sum(DIAGNOSIS == i), MinNARate = min(MinNARate, sum(is.na(Intensity) & DIAGNOSIS == i)/sum(DIAGNOSIS == i)), MaxNAcount = max(MaxNAcount, sum(is.na(Intensity) & DIAGNOSIS == i)))
    }
  }
  
  # then, filter based on percentage or raw amount
  if (method == "percentage") {
    combined_df <- combined_df %>%
      filter(MinNARate <= percentageNA)
  } else if (method == "raw") {
    combined_df <- combined_df %>%
      group_by(`Protein IDs`, DIAGNOSIS) %>%
      filter(MaxNAcount <= maxNA)
  } else {
    warning("Matrix was not filtered; method was not specified as 'percentage' or 'raw'")
  }
  
  return(combined_df)
}

# read protein group files
pgdata <- read_tsv("pgdata.tsv")
manifest <- read_tsv("manifest.tsv")
pgdata_processed <- runMQ(pgdata, log2t = TRUE)
# View(pgdata_processed)
# write_tsv(pgdata_processed, "pgdata_processed.tsv")
pgdata_2 <- matrix_filter(pgdata_processed, manifest, log2t = TRUE, method="raw", maxNA = 5)