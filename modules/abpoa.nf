process abpoa {
    tag { meta.well }
    label 'process_medium'
    publishDir "${params.out_dir}/consensus", mode: 'copy'

    conda (params.enable_conda ? 'bioconda::abpoa=1.5.4' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/abpoa:1.5.4--62e7dd62ad91bc78' :
        'community.wave.seqera.io/library/abpoa:1.5.4--cf15dcf6248e1556' }"

    input:
    tuple val(meta), path(heavy), path(light)

    output:
    tuple val(meta), path("*_consensus_heavy.fasta"), path("*_consensus_heavy.fasta"), path("*_consensus.fasta"), emit: consensus_seq

    script:

    """
    abpoa $heavy > "${meta.well}_consensus_heavy.fasta"
    sed -i '1c\\>${meta.well}_H' "${meta.well}_consensus_heavy.fasta"

    abpoa $light > "${meta.well}_consensus_light.fasta"
    sed -i '1c\\>${meta.well}_L' "${meta.well}_consensus_light.fasta"

    cat "${meta.well}_consensus_heavy.fasta" \
        "${meta.well}_consensus_light.fasta" > "${meta.well}_consensus.fasta"
    """
}