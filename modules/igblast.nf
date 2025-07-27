// module for running IgBLAST
// general layout is based on the nf-core modules
process igblast {
    tag { meta.well }
    label 'process_high'
    publishDir "${params.out_dir}/original_igblast", mode: 'copy'
        conda (params.enable_conda ? 'bioconda::igblast=1.19.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/igblast:7cbf716402092613' :
        'community.wave.seqera.io/library/igblast:9e55bf6484824037' }"

    input:
    tuple val(meta), path(reads), path(igblast_databases)
    val pre_post

    output:
    tuple val(meta), path('*_consensus_igblast.tsv'), emit: airr_table

    script:

    """
    # setup env variables
    export IGDATA="\${PWD}/igblast/igdata/"
    export IGBLASTDB="\${PWD}/igblast/databases"
    # run igblast
    # outfmt 19 = AIRR format (tsv, easy to use in downstream steps)
    igblastn -germline_db_V "\${PWD}/igblast/databases/imgt_human_V" \
        -germline_db_J "\${PWD}/igblast/databases/imgt_human_J" \
        -germline_db_D "\${PWD}/igblast/databases/imgt_human_D" \
        -organism human \
        -query $reads \
        -num_threads $task.cpus \
        -auxiliary_data "\${PWD}/igblast/igdata/optional_file/human_gl.aux" \
        -show_translation \
        -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 \
        -outfmt 19 > ${meta.well}_${pre_post}_consensus_igblast.tsv
    """
}