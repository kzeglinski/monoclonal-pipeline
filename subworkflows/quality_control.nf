// convert fastq to fasta (needed for IgBLAST)
process seqkit_stats {
    tag "qc_stats"
    label 'process_low'
    publishDir "${params.out_dir}/qc", mode: 'copy'

    conda (params.enable_conda ? 'bioconda::seqkit=2.3.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/seqkit:1c7b907eba99f587' :
        'community.wave.seqera.io/library/seqkit:825acf14813a21d5' }"

    input:
    path(fastq_files)

    output:
    path("*.tsv")

    script:
    """
    seqkit stats -a --threads $task.cpus ./* > seqkit_stats.tsv
    """
}

// adapted from the nf-core module: https://github.com/nf-core/modules/tree/master/modules/nf-core/minimap2/align
process minimap2_alignment {
    tag "$well"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::minimap2=2.30' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/minimap2:2.30--3bf3d6cb39a98dae' :
        'community.wave.seqera.io/library/minimap2:2.30--dde6b0c5fbc82ebd' }"

    input:
    tuple val(well), val(note), path(reads)
    path reference

    output:
    tuple val(well), path("*.sam"), emit: alignments

    script:
    """
    minimap2 \\
        -t $task.cpus \\
        -ax map-ont \\
        $reference \\
        $reads \\
        -o ${well}.sam
    """
}

process alignment_qc{
    tag "$well"
    label 'process_low'
    publishDir "${params.out_dir}/qc", mode: 'copy'

    conda (params.enable_conda ? 'bioconda::samtools=1.22.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/samtools:1.22.1--9a10f06c24cdf05f' :
        'community.wave.seqera.io/library/samtools:1.22.1--eccb42ff8fb55509' }"

    input:
    tuple val(well), path(alignments)

    output:
    tuple val(well), path("*.tsv"), emit: alignment_qc

    script:
    """
    # first need to sort 
    samtools view -b $alignments | samtools sort -o ${well}_sorted.bam

    # calculate the percent on target
    samtools flagstat -O tsv ${well}_sorted.bam > "${well}_flagstat.tsv"

    # then coverage
    samtools coverage ${well}_sorted.bam > "${well}_coverage.tsv"
 
    """
}

workflow irrelevant_qc{
    take:
        samples
        irrelevant_reference

    main:
        // specific qc for the irrelevant phage: align it to
        // the reference and check the coverage
        samples.map{it -> 
            def note = it[0].notes
            def well = it[0].well
            def reads = it[1]
            return [well: well, note: note, reads: reads]
        }
        .filter(it -> it.note == "irrelevant phage")
        .set{ irr_phage_samples }
        
        // mm2 alignment
        alignments = minimap2_alignment(irr_phage_samples, irrelevant_reference)

        // use samtools to get some qc stats on the alignment
        irr_phage_qc = alignment_qc(alignments)
    
    emit:
        irr_phage_samples
}

process post_consensus_qc {
    label 'process_low'
    publishDir "${params.out_dir}/key_outputs", mode: 'copy'

    conda (params.enable_conda ? 'conda-forge::r-tidyverse=2.0.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-biostrings_r-tidyverse:a28f97e8be230fb0' :
        'community.wave.seqera.io/library/bioconductor-biostrings_r-tidyverse:285d583c318ea12d' }"

    input:
    path(annotation_output)

    output:
    path("*.csv"), emit: consensus_qc
    path("*.fasta"), emit: final_fastas

    script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)

annotation <- read_tsv(

    # read in files from current task directory (contains files in annotation_output) and add column 'well' with the filename
    fs::dir_ls(glob = "*.tsv"), id = "well") %>%
    rename_with(~ "sequence_id", matches("sequence_(id|header)")) %>%
    select(sequence_id, locus, productive, v_call, 
        d_call, j_call, cdr1_aa, cdr2_aa, cdr3_aa, sequence, well) %>%
    
    # add count column (extract from sequence ID)
    mutate(count = as.numeric(str_remove(str_extract(sequence_id, "count_[[:digit:]]*"), "count_"))) %>% 
    relocate(count, .after = "sequence_id") %>%

    # add vector type column (esp. antibody vs nanobody) by extracting from the filename (well column)
    mutate(vector = case_when(
        str_detect(well, "nanobody") ~ "nanobody",
        str_detect(well, "antibody") ~ "antibody",
        str_detect(well, "empty_well") ~ "empty_well",
        str_detect(well, "irrelevant_phage") ~ "irrelevant_phage")) %>%
    
    # extract well from filename (remove everything after the first underscore which is after the well ID)
    mutate(well = str_remove(well, "_(?<=_).*")) %>%
    
    # for writing out the fasta later
    # mutate(sequence_id = paste0(">", str_remove(sequence_id, "_(?<=_).*"))) %>%
    
    # convert all values in 'locus' column to uppercase (riot outputs lowercase)
    mutate(locus = toupper(locus)) %>%

    # collapse those from the same well with same v&j
    # always take the productive one with highest count
    group_by(well, locus, v_call, cdr3_aa) %>%
    arrange(desc(productive), desc(count)) %>%
    summarise(
        well = well[1], 
        vector = vector[1],
        count = sum((count)), 
        locus = locus[1],
        productive = productive[1],
        v_call = v_call[1],
        d_call = d_call[1],
        j_call = j_call[1],
        cdr1_aa = cdr1_aa[1],
        cdr2_aa = cdr2_aa[1],
        cdr3_aa = cdr3_aa[1],
        sequence = sequence[1],
        sequence_id = sequence_id[1]
    )

fix_amber_stop <- function(annotation){
    amber_edit_code <- Biostrings::GENETIC_CODE
    amber_edit_code["TAG"] <- "Z" 

    # make dna string set from the sequence column
    dss <- annotation %>%
        pull(sequence) %>%
        Biostrings::DNAStringSet()

    # compare translation with regular + default code
    default_trans <- Biostrings::translate(dss)
    amber_trans <- Biostrings::translate(dss, genetic.code = amber_edit_code)

    code_comparison <- default_trans == amber_trans

    annotation <- annotation %>%
        ungroup() %>%
        mutate(
            has_amber_stop = !code_comparison,
            amber_stop_position = NA,
            amber_corrected_nt = NA
        )
    
    # loop over each row, cause idk how to use biostrings lol
    for (i in seq_len(nrow(annotation))) {
        # only need to do something if there was a diff 
        # between the default & amber translation
        if (!(code_comparison[i])) {
            # if amber stop(s) exist:
            # step 1: find them
            # this is the position of the T in TAG (* 3 to convert back
            # to nt and - 2 to get 1st nt in the codon 
            pos <- (BiocGenerics::start(Biostrings::matchPattern("Z", amber_trans[[i]])) * 3) - 2

            # replace it
            corrected_nt <- Biostrings::replaceAt(
                x = dss[[i]], # the untranslated DNAstring
                # need an IRanges of width 1 not just a numeric position otherwise
                # the replaceAt function does an insertion instead of substitution
                at = IRanges::IRanges(start = pos, width = 1),
                value = "C" # replace T with C
            ) 

            # now add this information to the annotation data 
            annotation[["amber_stop_position"]][i] <- pos
            annotation[["amber_corrected_nt"]][i] <- as.character(corrected_nt)
        }
    }
    return(annotation)
}

annotation <- fix_amber_stop(annotation)

# prepare output
# first thing is monoclonal qc table
monoclonal_qc <- data.frame(
    well = unique(annotation[["well"]]),
    vector = NA,
    monoclonal_status = NA
)

for (i in seq_along(unique(annotation[["well"]]))) {
    # get just this well's data
    this_well <- unique(annotation[["well"]])[i]
    this_annotation <- annotation %>%
        filter(well == this_well)

    # check for monoclonality

    # if vector type is empty_well or irrelevant_phage, skip monoclonal qc and just label as such
    
    if (this_annotation[["vector"]][1] == "empty_well") {
        status <- "empty_well"
    
    } else if (this_annotation[["vector"]][1] == "irrelevant_phage") {
        status <- "irrelevant_phage"
    
    # if vector type is antibody:
    
    } else if (this_annotation[["vector"]][1] == "antibody") {

        # if only 2 rows, it's monoclonal
        if (nrow(this_annotation) == 2) {
            status <- "monoclonal (no subclones)"
        
        # otherwise find counts of top/second H/L
        } else {

            top_heavy_count <- this_annotation %>% 
                filter(locus == "IGH") %>% arrange(desc(count)) %>% 
                dplyr::slice(1:1) %>% pull(count)

            if (nrow(filter(this_annotation, locus == "IGH")) > 1) {
                second_heavy_count <- this_annotation %>% 
                    filter(locus == "IGH") %>% arrange(desc(count)) %>% 
                    dplyr::slice(2:2) %>% pull(count)
            } else {
                second_heavy_count <- 0
            }

            top_light_count <- this_annotation %>% 
                filter(locus != "IGH") %>% arrange(desc(count)) %>% 
                dplyr::slice(1:1) %>% pull(count)

            if (nrow(filter(this_annotation, locus != "IGH")) > 1) {
                second_light_count <- this_annotation %>% 
                    filter(locus != "IGH") %>% arrange(desc(count)) %>% 
                    dplyr::slice(2:2) %>% pull(count)
            } else {
                second_light_count <- 0
            }
            
            # other chance to be monoclonal is if the second most abundant H/L 
            # chain is <10% the counts of the most abundant one

            if ((second_heavy_count <= 0.1 * top_heavy_count) & 
                (second_light_count <= 0.1 * top_light_count)) {
                    status <- "monoclonal (<10% subclones)"
                } else if ((second_heavy_count <= 0.2 * top_heavy_count) & 
                (second_light_count <= 0.2 * top_light_count)) {
                    status <- "likely monoclonal (<20% subclones)"
                } else if (top_heavy_count < 15 | top_light_count < 15) {
                    status <- "unsure (H or L count <15)"
                } else {
                    status <- "not monoclonal"
                }
        }
    
    # if vector type is nanobody:

    } else if (this_annotation[["vector"]][1] == "nanobody") {

        # if only 1 row, it's monoclonal
        if (nrow(this_annotation) == 1) {
            status <- "monoclonal"
        
        # otherwise find counts of top/second IGH clone
        # still filter for IGH just in case there are any non-IGH calls in the nanobody wells
        } else {
            top_heavy_count <- this_annotation %>%
                filter(locus == "IGH") %>% arrange(desc(count)) %>%
                dplyr::slice(1:1) %>% pull(count)

            if (nrow(this_annotation) > 1) {
                second_heavy_count <- this_annotation %>%
                    filter(locus == "IGH") %>% arrange(desc(count)) %>%
                    dplyr::slice(2:2) %>% pull(count)
            } else {
                second_heavy_count <- 0
            }

            if (second_heavy_count <= 0.1 * top_heavy_count) {
                status <- "monoclonal (<10% subclones)"
            } else if (second_heavy_count <= 0.2 * top_heavy_count) {
                status <- "likely monoclonal (<20% subclones)"
            } else if (top_heavy_count < 15) {
                status <- "unsure (H count <15)"
            } else {
                status <- "not monoclonal"
            }
        }
    }

    monoclonal_qc[["vector"]][i] <- this_annotation[["vector"]][1]
    monoclonal_qc[["monoclonal_status"]][i] <- status

}

# now, write things out
# monoclonal table
write_csv(monoclonal_qc, "monoclonal_qc.csv")

# annotation table
# remove vector column from monoclonal_qc for joining
monoclonal_qc <- monoclonal_qc %>% select(-vector)
annotation %>% 
    left_join(monoclonal_qc, by = "well") %>%
    write_csv("consensus_annotation.csv")

# fasta with all subclones included
all_subclone_fasta <- annotation %>%
    # create read headers
    mutate(chain = case_when(locus == "IGH" ~ "H", TRUE ~ "L")) %>%
    mutate(read_name = paste0(
        ">", well, "_", vector, "_", chain, "_count:", count
    )) %>%
    # combine read name and seq
    mutate(fasta = paste(read_name, sequence, sep = "\n")) %>%
    pull(fasta)

write_lines(all_subclone_fasta, "all_subclones.fasta")

# fasta with amber correction and no subclones
amber_corrected_fasta <- annotation %>%
    # grab only top clones
    group_by(well, locus) %>%
    slice_max(count) %>%
    ungroup() %>%
    # create read headers
    mutate(
        chain = case_when(locus == "IGH" ~ "H", TRUE ~ "L"),
        corrected = ifelse(has_amber_stop, "_amber_corrected", ""),
        sequence = case_when(
            !is.na(amber_corrected_nt) ~ amber_corrected_nt,
            TRUE ~ sequence)) %>%
    mutate(read_name = paste0(
        ">", well, "_", vector, "_", chain, "_count:", count, corrected
    )) %>%
    # combine read name and seq
    mutate(fasta = paste(read_name, sequence, sep = "\n")) %>%
    pull(fasta)

write_lines(amber_corrected_fasta, "top_clones.fasta")

    """
}