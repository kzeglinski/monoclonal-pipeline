include { igblast } from '../modules/igblast'

process group_reads {
    tag { meta.well }
    label 'process_low'
    publishDir "${params.out_dir}/qc", mode: 'copy', pattern: "*.tsv"

    conda (params.enable_conda ? 'conda-forge::r-tidyverse=2.0.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-tidyverse:2.0.0--33f4d800f6070aac' :
        'community.wave.seqera.io/library/r-tidyverse:2.0.0--dd61b4cbf9e28186' }"

    input:
    tuple val(meta), path(igblast_output)

    output:
    tuple val(meta), path("*_heavy_clean.fasta"), path("*_light_clean.fasta"), emit: consensus_input, optional: true
    path("*_qc.tsv"), emit: qc

    script:
    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    id <- '$meta.well'
    annotation <- read_tsv(
        fs::dir_ls(glob = "*_pre_consensus_igblast.tsv"),
        col_select = c(sequence_id, complete_vdj, v_call, 
            d_call, j_call, sequence)) %>%
        # for writing out the fasta later
        mutate(sequence_id = paste0(">", sequence_id)) %>%
        # remove abnormally long sequences, they mess up consensus
        filter(nchar(sequence) < 500)

    # heavy chains
    heavy <- annotation %>%
        filter(str_detect(v_call, "IGH")) %>%
        group_by(v_call, d_call, j_call) %>%
        summarise(
            count = n(),
            reads = paste(sequence_id, sequence, 
                sep = "\n", collapse = "\n")) %>%
        arrange(desc(count))

    light <- annotation %>%
        filter(!str_detect(v_call, "IGH")) %>%
        group_by(v_call, d_call, j_call) %>%
        summarise(
            count = n(),
            reads = paste(sequence_id, sequence, 
                sep = "\n", collapse = "\n")) %>%
        arrange(desc(count))

    # check if monoclonal -> 2nd most abudant v/d/j combo should
    # be 10% or less of the main one
    main_heavy <- ungroup(heavy) %>% slice_head(n = 1) %>% pull(count)
    main_light <- ungroup(light) %>% slice_head(n = 1) %>% pull(count)

    if (nrow(light) == 1){
        second_light <- 0
    } else {
        second_light <- ungroup(light) %>% slice(2:2) %>% pull(count)
    }

    if (nrow(heavy) == 1){
        second_heavy <- 0
    } else {
        second_heavy <- ungroup(heavy) %>% slice(2:2) %>% pull(count)
    }

    if ((main_heavy * 0.15 < second_heavy) | (main_light * 0.15 < second_light)){
        status <- "not monoclonal"
    } else {
        status <- "normal"
        # write out the heavy and light reads for consensus
        ungroup(heavy) %>% 
            slice_head(n = 1) %>% 
            pull(reads) %>%
            write_lines(file = paste0(id, "_heavy_clean.fasta"))

        ungroup(light) %>% 
            slice_head(n = 1) %>% 
            pull(reads) %>%
            write_lines(file = paste0(id, "_light_clean.fasta"))
    }

    monoclonal_status <- data.frame(
        id = id,
        status = status,
        heavy_counts = paste(main_heavy, second_heavy, sep = "_"),
        light_counts = paste(main_light, second_light, sep = "_")
    )
    write_tsv(monoclonal_status, paste0(id, "_monoclonal_qc.tsv"))
    """
}

workflow pre_consensus_grouping {
    take:
        reads_with_igblast_data

    main:
        // annotate reads using igblast
        // TO-DO: GET USER TO INPUT ORGANISM SOMEWHERE? IN FLANK FILE?
        // MAYBE WE COULD RENAME FLANK FILE LIKE LIB STRUCT OR SMTH
        // ATM ITS HARDCODED TO HUMAN

        igblast_tsv = igblast(
            reads_with_igblast_data,
            "pre").airr_table

        // process this in R
        grouped_reads = group_reads(igblast_tsv)
        consensus_input = grouped_reads.consensus_input
        qc = grouped_reads.qc
    emit:
        consensus_input
        qc
}