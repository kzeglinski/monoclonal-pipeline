process preconsensus_group_reads {
    tag { meta.well }
    label 'process_low'
    publishDir "${params.out_dir}/antibody_clone_counts", pattern: "*_antibody_clones.csv", mode: 'copy'

    conda (params.enable_conda ? 'conda-forge::r-tidyverse=2.0.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-tidyverse:2.0.0--33f4d800f6070aac' :
        'community.wave.seqera.io/library/r-tidyverse:2.0.0--dd61b4cbf9e28186' }"

    input:
    tuple val(meta), path(annotation_output)

    output:
    tuple val(meta), path("*_clean.fasta"), emit: consensus_input, optional: true
    path("*_antibody_clones.csv"), emit: antibody_clone_counts, optional: true

    script:
    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    id <- '$meta.well'
    vector <- '$meta.vector'
    
    annotation <- read_tsv(
        fs::dir_ls(glob = "*pre_consensus_annotation.tsv")) %>%
        rename_with(~ "sequence_id", matches("sequence_(id|header)")) %>%
        select(sequence_id, complete_vdj, v_call, d_call, j_call, sequence) %>%
        # for writing out the fasta later
        mutate(sequence_id = paste0(">", sequence_id)) %>%
        # remove abnormally long sequences, they mess up consensus
        filter(nchar(sequence) < 1000)
    
    # only proceed if there are any annotated reads (i.e. has v_call, d_call, and j_call)
    if (nrow(annotation %>% filter(!is.na(v_call) & !is.na(d_call) & !is.na(j_call))) > 0) {
        heavy <- annotation %>%
            filter(str_detect(v_call, "IGH|VHH|vhh")) %>%
            summarise(
            # count and sequences (fasta format)
            # for reads with the same IGH v/d/j call combo (i.e. belong to the same clone)
            # arranged in descending order by count (i.e. most abundant clone first)
                count = n(),
                reads = paste(sequence_id, sequence, 
                    sep = "\n", collapse = "\n"),
                .by = c(v_call, d_call, j_call)) %>%
            arrange(desc(count))

        # new approach: take consensus of all v/d/j combo that is has a count 
        # at least 10% of the most abundant combo
        ## change to 5% to provide some leeway for heavy/light chain clone combos
        ## may change in future to just take into account heavy/light chain clone pairs

        # first find that cutoff
        # cons_cutoff_h <- floor(pull(slice_head(heavy, n = 1), count) * 0.1)
        cons_cutoff_h <- floor(pull(slice_head(heavy, n = 1), count) * 0.05)

        # count needs to be at least 3
        # if not, just take the most abundant one 
        if (cons_cutoff_h < 3) {
            cons_cutoff_h <- pull(slice_head(heavy, n = 1), count)
        }

        # write out the heavy reads for consensus
        # do it for any that pass the cutoff
        for_cons_h <- heavy %>% filter(count >= cons_cutoff_h)

        for (i in seq_len(nrow(for_cons_h))) {
            this_clone <- for_cons_h[i ,]
            this_count <- this_clone\$count
            this_clone %>% pull(reads) %>% write_lines(file = paste0(id, "_", i , "_count_", this_count, "_heavy_clean.fasta"))
        }

        # if vector is not nanobody (i.e. it's antibody or empty_well or irrelevant_phage), also write out the light chains for consensus

        if (vector != "nanobody") {
            light <- annotation %>%
                    filter(!str_detect(v_call, "IGH")) %>%
                    summarise(
                        count = n(),
                        reads = paste(sequence_id, sequence, 
                            sep = "\n", collapse = "\n"),
                        .by = c(v_call, j_call)) %>%
                    arrange(desc(count))
            
            # take consensus of all v/d/j combo that is has a count 
            # at least 10% of the most abundant combo
            ## change to 5% to provide some leeway for heavy/light chain clone combos
            ## may change in future to just take into account heavy/light chain clone pairs

            # first find that cutoff
            # cons_cutoff_l <- floor(pull(slice_head(light, n = 1), count) * 0.1)
            cons_cutoff_l <- floor(pull(slice_head(light, n = 1), count) * 0.05)
            
            # count needs to be at least 3
            # if not, just take the most abundant one 

            if (cons_cutoff_l < 3) {
                cons_cutoff_l <- pull(slice_head(light, n = 1), count)
            }

            # write out the light reads for consensus
            # do it for any that pass the cutoff
            for_cons_l <- light %>% filter(count >= cons_cutoff_l)
            
            for (i in seq_len(nrow(for_cons_l))) {
                this_clone <- for_cons_l[i ,]
                this_count <- this_clone\$count
                this_clone %>% pull(reads) %>% write_lines(file = paste0(id, "_", i , "_count_", this_count,"_light_clean.fasta"))
            }

            # get unique heavy v/d/j call combos and assign clone name (v_call + d_call + j_call)
            heavy_clones <- annotation %>%
                filter(str_detect(v_call, "IGH|VHH|vhh")) %>%
                select(v_call, d_call, j_call) %>%
                distinct() %>%
                mutate(heavy_clone = paste0(v_call, "_", d_call, "_", j_call))
            
            # get unique light v/j call combos and assign clone name (v_call + j_call)
            light_clones <- annotation %>%
                filter(!str_detect(v_call, "IGH")) %>%
                select(v_call, j_call) %>%
                distinct() %>%
                mutate(light_clone = paste0(v_call, "_", j_call))
            
            # annotate each sequence_id with its heavy clone and light clone based on its v/d/j calls
            clone_annotation <- annotation %>%
                left_join(heavy_clones, by = c("v_call", "d_call", "j_call")) %>%
                left_join(light_clones, by = c("v_call", "j_call")) %>%
                select(sequence_id, v_call, d_call, j_call, heavy_clone, light_clone)
            
            # combine rows with same sequence_id (i.e. same read) - summarise .by(sequence_id)
            # tibble will have rows:
            # sequence_id, heavy_clone, light_clone

            read_clone_annotation <- clone_annotation %>%
                group_by(sequence_id) %>%
                summarise(
                    # only include heavy_clone and light_clone that are not NA
                    heavy_clone = paste(unique(heavy_clone[!is.na(heavy_clone)]), collapse = ";"),
                    light_clone = paste(unique(light_clone[!is.na(light_clone)]), collapse = ";")
                )
            
            # now count how many reads belong to each unique combo of heavy_clone and light_clone
            clone_combos <- read_clone_annotation %>%
                group_by(heavy_clone, light_clone) %>%
                summarise(count = n()) %>%
                arrange(desc(count)) %>%
                ungroup()
            
            # print out the clone combos and their counts to a csv
            write_csv(clone_combos, paste0(id, "_antibody_clones.csv"))
        }
    }
    
    """
}
