# install packages
library(tidyverse)
#install.packages("fuzzyjoin")
library(fuzzyjoin)
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("pcaMethods")
library(pcaMethods)
# library(imputeLCMD)

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

combine_matrix <- function (pg_df, mf_df, sample.id.col='Sample.ID', intensity_type = 'Raw') {
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
    mutate(Sample.ID.regex = paste("(Intensity|LFQ intensity|MS/MS count)\\s", !!sym(sample.id.col), ".*", sep=""))
  
  gsub(".", "[.]", manifest_local$Sample.ID.regex)
  
  # use left join to add patient diagnosis to sample ID
  combined_df <- regex_left_join(proteingroups, manifest_local, by=c(Sample = "Sample.ID.regex"))[-ncol(proteingroups) - ncol(manifest_local)] %>%
    select(-!!sym(sample.id.col)) %>%
    distinct()
  
  return(combined_df)
}

filter_intensities <- function(pg_df, method='percentage', filterGroup='DIAGNOSIS',percentageNA = 0.5, maxNA=0, log2t) {
  # change -Inf to NA
  is.na(pg_df) <- sapply(pg_df, is.infinite)
  
  # throw error if percentageNA < 0 or > 100
  if (percentageNA < 0 || percentageNA > 100) {
    stop("percentageNA must be between 0 and 100.")
  }
  
  # create NARate column names for every unique diagnosis
  diagnoses <- unique(pull(pg_df, !!sym(filterGroup)))
  
  # filter values that do not meet conditions
  # first, create columns of NA counts and total diagnoses for each protein group and diagnosis
  pg_df <- pg_df %>%
    mutate(MinNARate = 1, MaxNAcount = 0)
  
  if (log2t) {
    for (i in diagnoses) {
      pg_df <- pg_df %>%
        group_by(`Protein IDs`) %>%
        mutate("NACount.{{i}}" := sum(is.na(Intensity) & !!sym(filterGroup) == i), "NARate.{{i}}" := sum(is.na(Intensity) & !!sym(filterGroup) == i)/sum(!!sym(filterGroup) == i), MinNARate = min(MinNARate, sum(is.na(Intensity) & !!sym(filterGroup) == i)/sum(!!sym(filterGroup) == i)), MaxNAcount = max(MaxNAcount, sum(is.na(Intensity) & !!sym(filterGroup) == i)))
    }
  } else {
    for (i in diagnoses) {
      pg_df <- pg_df %>%
        group_by(`Protein IDs`) %>%
        mutate("NACount.{{i}}" := sum(Intensity == 0 & !!sym(filterGroup) == i), "NARate.{{i}}" := sum(is.na(Intensity) & !!sym(filterGroup) == i)/sum(!!sym(filterGroup) == i), MinNARate = min(MinNARate, sum(is.na(Intensity) & !!sym(filterGroup) == i)/sum(!!sym(filterGroup) == i)), MaxNAcount = max(MaxNAcount, sum(is.na(Intensity) & !!sym(filterGroup) == i)))
    }
  }
  
  # then, filter based on percentage or raw amount
  if (method == "percentage") {
    pg_df <- pg_df %>%
      filter(MinNARate <= percentageNA)
  } else if (method == "raw") {
    pg_df <- pg_df %>%
      group_by(`Protein IDs`, !!sym(filterGroup)) %>%
      filter(MaxNAcount <= maxNA)
  } else {
    warning("Matrix was not filtered; method was not specified as 'percentage' or 'raw'")
  }
  
  return(pg_df)
}

impute <- function(pg_df, method = "PPCA", nPcs=3, imputeGroup = "DIAGNOSIS", intensity_type) {
  # tibble --> data matrix
  
  # split into groups
  diagnoses <- unique(pull(pg_df, !!sym(imputeGroup)))
  
  groupdata_df <- pg_df %>%
    select(`Protein IDs`, `Gene names`, Reverse, `Potential contaminant`, id, `Peptide IDs`, Sample, Intensity)
  impute_df <- tibble(groupdata_df) %>%
    ungroup() %>%
    select(-!!sym(imputeGroup)) %>%
    spread(Sample, Intensity) %>%
    select(-starts_with("Intensity"), -starts_with("LFQ intensity"), -starts_with("MS/MS count"))
  # impute each group, combine with impute_df
  for (i in diagnoses) {
    cur_grp <- groupdata_df %>%
      ungroup() %>%
      filter(!!sym(imputeGroup) == i) %>%
      spread(Sample, Intensity) %>%
      select(`Protein IDs`, starts_with("Intensity"), starts_with("LFQ intensity"), starts_with("MS/MS count"), id)
    pc <- pca(cur_grp, method="bpca", nPcs=nPcs)
    impute_df <- impute_df %>%
      inner_join(completeObs(pc), copy=TRUE)
  }
  
  # convert to long format
  if (intensity_type == "Raw") {
    impute_df <- impute_df %>%
      gather("Sample", "Intensity", starts_with("Intensity"), starts_with("LFQ intensity"), starts_with("MS/MS count"))
  } else if (intensity_type == "LFQ") {
    impute_df <- impute_df %>%
      gather("Sample", "LFQ intensity", starts_with("Intensity"), starts_with("LFQ intensity"), starts_with("MS/MS count"))
  } else if (intensity_type == "MS2") {
    impute_df <- impute_df %>%
      gather("Sample", "MS/MS count", starts_with("Intensity"), starts_with("LFQ intensity"), starts_with("MS/MS count"))
  } else {
    warning("type should be specified as 'Raw', 'LFQ', or 'MS2'.")
  }
  
  return (impute_df)
}

# read protein group files
pgdata <- read_tsv("pgdata.tsv")
manifest <- read_tsv("manifest.tsv")
pgdata_1 <- runMQ(pgdata, log2t = TRUE)
pgdata_2 <- combine_matrix(pgdata_processed, manifest, intensity_type="Raw")
pgdata_3 <- filter_intensities(pgdata_2, log2t = TRUE, method="raw", maxNA = 5)
pgdata_4 <- impute(pgdata_3, nPcs = 3, intensity_type="Raw")
write_tsv(pgdata_4, "pgdata_processed.tsv")