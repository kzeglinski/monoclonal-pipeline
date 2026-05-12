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

# check number of flanking sequences for each vector type
nanobody_nb_rows <- nrow(flanking_sequences[flanking_sequences[["vector_type"]] == "nanobody" & 
                                            flanking_sequences[["flank_name"]] == "nb", ])
antibody_vh_rows <- nrow(flanking_sequences[flanking_sequences[["vector_type"]] == "antibody" & 
                                            flanking_sequences[["flank_name"]] == "vh", ])
antibody_vlk_rows <- nrow(flanking_sequences[flanking_sequences[["vector_type"]] == "antibody" & 
                                            flanking_sequences[["flank_name"]] == "vlk", ])
antibody_vll_rows <- nrow(flanking_sequences[flanking_sequences[["vector_type"]] == "antibody" &
                                            flanking_sequences[["flank_name"]] == "vll", ])

# if vector_type is 'nanobody', check that there is only 1 nanobody nb flanking sequence
if(vector_type == "nanobody"){
    if(nanobody_nb_rows != 1){
        stop("For 'nanobody' vector type, the flanking sequences file should contain 1 'nb' flanking sequence. Found ", 
            nanobody_nb_rows, " 'nb' nanobody flanking sequences.")
    }
}

# if vector_type is 'antibody', check that there is only 1 antibody vh, vlk and vll flanking sequence
if(vector_type == "antibody"){
    if(antibody_vh_rows != 1 || antibody_vlk_rows != 1 || antibody_vll_rows != 1){
        stop("For 'antibody' vector type, the flanking sequences file should contain 1 'vh' sequence, 1 'vlk' sequence and 1 'vll' sequence.\n Found ",
            antibody_vh_rows, " 'vh' antibody sequence(s), ",
            antibody_vlk_rows, " 'vlk' antibody sequence(s) and ",
            antibody_vll_rows, " 'vll' antibody sequence(s).")
    }
}

# if vector_type is 'both', check that there is only 1 nanobody nb flanking sequence and 1 antibody vh, vlk and vll flanking sequence
if(vector_type == "both"){
    if(nanobody_nb_rows != 1 || antibody_vh_rows != 1 || antibody_vlk_rows != 1 || antibody_vll_rows != 1){
        stop("For 'both' vector type (antibodies & nanobodies), the flanking sequences file should contain 1 'nb' nanobody sequence, 1 'vh' antibody sequence, 1 'vlk' antibody sequence and 1 'vll' antibody sequence.\n Found ",
            nanobody_nb_rows, " 'nb' nanobody sequence(s), ",
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

flanks <- list()
for(i in 1:nrow(flanking_sequences)){
    row <- flanking_sequences[i, ]
    # join vector type, flank name and L/R to create parameter name (e.g. vh_flank_L, vlk_flank_R, etc.)
    flank_name <- paste0(row[["vector_type"]], "_", row[["flank_name"]])
    flanks[[paste0(flank_name, "_L")]] <- row[["flank_L"]]
    flanks[[paste0(flank_name, "_R")]] <- row[["flank_R"]]
}

flanks_df <- data.frame(
    parameter = names(flanks),
    sequence = unlist(flanks)
)

write_csv(flanks_df, "flanking_seqs_for_mb.csv")