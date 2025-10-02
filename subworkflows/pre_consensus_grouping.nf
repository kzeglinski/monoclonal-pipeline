include { igblast } from '../modules/igblast'

process group_reads {
    tag { meta.well }
    label 'process_low'

    conda (params.enable_conda ? 'conda-forge::r-tidyverse=2.0.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-tidyverse:2.0.0--33f4d800f6070aac' :
        'community.wave.seqera.io/library/r-tidyverse:2.0.0--dd61b4cbf9e28186' }"

    input:
    tuple val(meta), path(igblast_output)

    output:
    tuple val(meta), path("*_heavy_clean.fasta"), path("*_light_clean.fasta"), emit: consensus_input, optional: true

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
        filter(nchar(sequence) < 1000)

    # heavy chains
    heavy <- annotation %>%
        filter(str_detect(v_call, "IGH")) %>%
        summarise(
            count = n(),
            reads = paste(sequence_id, sequence, 
                sep = "\n", collapse = "\n"),
            .by = c(v_call, d_call, j_call)) %>%
        arrange(desc(count))

    light <- annotation %>%
        filter(!str_detect(v_call, "IGH")) %>%
        summarise(
            count = n(),
            reads = paste(sequence_id, sequence, 
                sep = "\n", collapse = "\n"),
            .by = c(v_call, j_call)) %>%
        arrange(desc(count))

    # new approach: take consensus of all v/d/j combo that is has a count 
    # at least 10% of the most abundant combo

    # first find that cutoff
    cons_cutoff_h <- floor(pull(slice_head(heavy, n = 1), count) * 0.1)
    cons_cutoff_l <- floor(pull(slice_head(light, n = 1), count) * 0.1)

    # count needs to be at least 3
    # if not, just take the most abundant one 
    if (cons_cutoff_h < 3) {
        cons_cutoff_h <- pull(slice_head(heavy, n = 1), count)
    }

    if (cons_cutoff_l < 3) {
        cons_cutoff_l <- pull(slice_head(light, n = 1), count)
    }

    # write out the heavy and light reads for consensus
    # do it for any that pass the cutoff
    for_cons_h <- heavy %>% filter(count >= cons_cutoff_h)
    for_cons_l <- light %>% filter(count >= cons_cutoff_l)

    for (i in seq_len(nrow(for_cons_h))) {
        this_clone <- for_cons_h[i ,]
        this_count <- this_clone\$count
        this_clone %>% pull(reads) %>% write_lines(file = paste0(id, "_", i , "_count_", this_count, "_heavy_clean.fasta"))
    }

    for (i in seq_len(nrow(for_cons_l))) {
        this_clone <- for_cons_l[i ,]
        this_count <- this_clone\$count
        this_clone %>% pull(reads) %>% write_lines(file = paste0(id, "_", i , "_count_", this_count,"_light_clean.fasta"))
    }
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

    emit:
        consensus_input

}