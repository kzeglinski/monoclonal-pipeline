#!/usr/bin/env Rscript
# sample_sheet_validation.r

library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
sample_sheet_file <- args[1]

# check that sample sheet exists
if(!file.exists(sample_sheet_file)){
    stop("Sample sheet does not exist at ", sample_sheet_file, "\n Please check the path and try again.")
}

# read it in
sample_sheet <- read_csv(sample_sheet_file, show_col_types = FALSE)

# check all required columns present
required_columns <- c("barcode", "well", "sample_id", "vector")
missing_columns <- setdiff(required_columns, colnames(sample_sheet))

if(length(missing_columns) > 0){
    stop("Sample sheet is missing the following required column(s): ", 
        paste(missing_columns, collapse = ", "))
}

# check all barcodes are of the form 'barcode01, barcode02, ..., barcode96'
barcode_pattern <- "^barcode(0[1-9]|[1-8][0-9]|9[0-6])$"
invalid_barcodes <- sample_sheet[["barcode"]][!grepl(barcode_pattern, sample_sheet[["barcode"]])]

if(length(invalid_barcodes) > 0){
    stop("The 'barcode' column contains the following invalid barcode(s): ", 
        paste(invalid_barcodes, collapse = ", "),
        ". Barcodes should be between 'barcode01' and 'barcode96'.")
}

# check for duplicate barcodes
duplicated_barcodes <- sample_sheet[["barcode"]][duplicated(sample_sheet[["barcode"]])]

if(length(duplicated_barcodes) > 0){
    stop("The 'barcode' column contains duplicate barcode(s): ", 
        paste(duplicated_barcodes, collapse = ", "),
        ". Each barcode must be unique.")
}

# check all well values are the correct format (letter from A to H followed by one of 01, 02, ..., 12)
well_pattern <- "^[A-H](0[1-9]|1[0-2])$"
invalid_wells <- sample_sheet[["well"]][!grepl(well_pattern, sample_sheet[["well"]])]

if (length(invalid_wells) > 0){
    stop("The 'well' column contains the following invalid well(s): ", 
        paste(invalid_wells, collapse = ", "),
        ". Wells should be in the format of a letter from A to H followed by a number from 01 to 12 (e.g. A01, B01, ..., H12).")
}

# check for duplicate wells
duplicated_wells <- sample_sheet[["well"]][duplicated(sample_sheet[["well"]])]

if(length(duplicated_wells) > 0){
    stop("The 'well' column contains duplicate well(s): ", 
        paste(duplicated_wells, collapse = ", "),
        ". Each well must be unique.")
}

# check vector column
sample_sheet <- sample_sheet %>%
    mutate(vector = case_when(
        str_detect(vector, regex("empty", ignore_case = TRUE)) ~ "empty_well",
        str_detect(vector, regex("irr", ignore_case = TRUE)) ~ "irrelevant_phage",
        str_detect(vector, regex("anti", ignore_case = TRUE)) | str_detect(vector, regex("ab", ignore_case = TRUE)) ~ "antibody",
        str_detect(vector, regex("nano", ignore_case = TRUE)) | str_detect(vector, regex("nb", ignore_case = TRUE)) ~ "nanobody",
        TRUE ~ str_replace_all(tolower(vector), "\\s+", "")
    ))

allowed_vector_values <- c("antibody", "nanobody", "empty_well", "irrelevant_phage")
invalid_vector_values <- sample_sheet[["vector"]][!sample_sheet[["vector"]] %in% allowed_vector_values]

if(length(invalid_vector_values) > 0){
    stop("The 'vector' column contains the following invalid value(s): ", 
        paste(invalid_vector_values, collapse = ", "),
        ". Allowed values are: antibody, nanobody, empty_well, irrelevant_phage.")
}

write_csv(sample_sheet, "validated_sample_sheet.csv")