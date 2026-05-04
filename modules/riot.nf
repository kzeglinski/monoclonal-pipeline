// module for running riot

process riot {
    tag { meta.well }
    label 'process_high'
        conda (params.enable_conda ? 'bioconda::riot-na=4.0.2' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/biopython_gcc_libxcrypt_python_pruned:96826a8c3e510274' :
        'community.wave.seqera.io/library/biopython_gcc_libxcrypt_python_pruned:c0ab77d048c45319' }"
    
    input:
    tuple val(meta), path(reads)
    val pre_post

    output:
    tuple val(meta), path('*_consensus_annotation.csv'), emit: airr_table

    script:
    def species = meta.vector == 'nanobody' ? 'VICUGNA_PACOS' : 'HOMO_SAPIENS'

    """
    riot_na -f $reads --species ${species} -p 16 -o "${meta.well}_${meta.vector}_${pre_post}_consensus_annotation.csv"
    """
}

// convert to tsv to make consistent with igblast output
// for input into preconsensus_group_reads()

process csv_to_tsv {
    tag "csv_to_tsv"
    label 'process_low'
    publishDir "${params.out_dir}/original_riot", mode: 'copy', pattern: "*_consensus_annotation.tsv"

    conda (params.enable_conda ? 'conda-forge::r-tidyverse=2.0.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-tidyverse:2.0.0--33f4d800f6070aac' :
        'community.wave.seqera.io/library/r-tidyverse:2.0.0--dd61b4cbf9e28186' }"

    input:
    tuple val(meta), path(annotation_output)

    output:
    tuple val(meta), path("*_consensus_annotation.tsv"), emit: tsv

    script:
    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    
    annotation_csv <- read_csv(fs::dir_ls(glob = "*_consensus_annotation.csv"))
    annotation_tsv <- write_tsv(annotation_csv, paste0(tools::file_path_sans_ext(fs::dir_ls(glob = "*_consensus_annotation.csv")), ".tsv"))
    """
}