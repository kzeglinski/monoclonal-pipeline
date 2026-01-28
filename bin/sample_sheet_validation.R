#!/usr/bin/env Rscript
# sample_sheet_validation.r

library(readr)

args <- commandArgs(trailingOnly=TRUE)
sample_sheet_file <- args[1]

# check that sample sheet exists
if(!file.exists(sample_sheet_file)){
    stop("Sample sheet does not exist at ", sample_sheet_file, "\n Please check the path and try again.")
}

# read it in
sample_sheet <- read_csv(sample_sheet_file, show_col_types = FALSE)

# check all required columns present
required_columns <- c("barcode", "alias", "cbt_phage_id", "notes")
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

# check all alias values are the correct format (letter from A to H followed by one of 01, 02, ..., 12)
alias_pattern <- "^[A-H](0[1-9]|1[0-2])$"
invalid_aliases <- sample_sheet[["alias"]][!grepl(alias_pattern, sample_sheet[["alias"]])]

if (length(invalid_aliases) > 0){
    stop("The 'alias' column contains the following invalid alias(es): ", 
        paste(invalid_aliases, collapse = ", "),
        ". Aliases should be in the format of a letter from A to H followed by a number from 01 to 12 (e.g. A01, B01, ..., H12).")
}

# check for duplicate aliases
duplicated_aliases <- sample_sheet[["alias"]][duplicated(sample_sheet[["alias"]])]

if(length(duplicated_aliases) > 0){
    stop("The 'alias' column contains duplicate alias(es): ", 
        paste(duplicated_aliases, collapse = ", "),
        ". Each alias must be unique.")
}