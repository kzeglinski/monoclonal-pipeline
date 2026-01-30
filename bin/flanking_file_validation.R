#!/usr/bin/env Rscript
# flanking_file_validation.R

library(readr)

args <- commandArgs(trailingOnly=TRUE)
flanking_sequences_file <- args[1]
vector_type <- args[2]

# check that flanking file exists
if(!file.exists(flanking_sequences_file)){
    stop("Flanking sequences file does not exist at ", flanking_sequences_file, "\n Please check the path and try again.")
}

# read it in
flanking_sequences <- read_csv(flanking_sequences_file, show_col_types = FALSE)

# check all columns present
required_columns <- c("vector_type", "flank_name", "flank_L", "flank_R")
missing_columns <- setdiff(required_columns, colnames(flanking_sequences))

if(length(missing_columns) > 0){
    stop("Flanking sequences file is missing the following required column(s): ", 
        paste(missing_columns, collapse = ", "))
}

# if vector_type is 'nanobody', check that there is only one nanobody vh flanking sequence
if(vector_type == "nanobody"){
    nanobody_vh_rows <- nrow(flanking_sequences[flanking_sequences[["vector_type"]] == "nanobody" & 
                                                flanking_sequences[["flank_name"]] == "vh", ])
    if(nanobody_vh_rows != 1){
        stop("For 'nanobody' vector type, the flanking sequences file should contain 1 'vh' flanking sequence. Found ", 
            nanobody_vh_rows, " 'vh' nanobody flanking sequences.")
    }
}

# if vector_type is 'antibody', check that there is only one antibody vh, vlk and vll flanking sequence
if(vector_type == "antibody"){
    antibody_flank_names <- c("vh", "vlk", "vll")
    antibody_vh_rows <- nrow(flanking_sequences[flanking_sequences[["vector_type"]] == "antibody" & 
                                                flanking_sequences[["flank_name"]] == "vh", ])
    antibody_vlk_rows <- nrow(flanking_sequences[flanking_sequences[["vector_type"]] == "antibody" & 
                                                flanking_sequences[["flank_name"]] == "vlk", ])
    antibody_vll_rows <- nrow(flanking_sequences[flanking_sequences[["vector_type"]] == "antibody" &
                                                flanking_sequences[["flank_name"]] == "vll", ])
    if(antibody_vh_rows != 1 || antibody_vlk_rows != 1 || antibody_vll_rows != 1){
        stop("For 'antibody' vector type, the flanking sequences file should contain 1 'vh' sequence, 1 'vlk' sequence and 1 'vll' sequence.\n Found ",
            antibody_vh_rows, " 'vh' antibody sequence(s), ",
            antibody_vlk_rows, " 'vlk' antibody sequence(s) and ",
            antibody_vll_rows, " 'vll' antibody sequence(s).")
    }
}

# check that flank_L and flank_R contain only valid DNA characters (A, T, C, G)
dna_pattern <- "^[ATCGatcg]+$"
invalid_flank_L <- flanking_sequences[["flank_L"]][!grepl(dna_pattern, flanking_sequences[["flank_L"]])]
invalid_flank_R <- flanking_sequences[["flank_R"]][!grepl(dna_pattern, flanking_sequences[["flank_R"]])]

if(length(invalid_flank_L) > 0){
    stop("The 'flank_L' column contains the following invalid sequence(s): ", 
        paste(invalid_flank_L, collapse = ", "),
        ". Sequences should contain only DNA characters (A, T, C, G).")
}

if(length(invalid_flank_R) > 0){
    stop("The 'flank_R' column contains the following invalid sequence(s): ", 
        paste(invalid_flank_R, collapse = ", "),
        ". Sequences should contain only DNA characters (A, T, C, G).")
}