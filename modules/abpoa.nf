process abpoa {
    tag { meta.well }
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::abpoa=1.5.4' : null)
    // container with just abpoa
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'oras://community.wave.seqera.io/library/abpoa:1.5.4--62e7dd62ad91bc78' :
    //    'community.wave.seqera.io/library/abpoa:1.5.4--cf15dcf6248e1556' }"
    
    // container with abpoa and seqkit for subsetting reads first
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'oras://community.wave.seqera.io/library/abpoa_seqkit:9b947cb71c5cc140' :
    'community.wave.seqera.io/library/abpoa_seqkit:0e9bcb3f2f1c3292' }"

    input:
    tuple val(meta), path(heavy), path(light)

    output:
    tuple val(meta), path("*_all_consensus.fasta"), emit: consensus_seq

    script:

    """
    for file in *_heavy_clean.fasta; do
        # Skip if no files matched
        [ -e "\$file" ] || continue

        base="\$(basename "\$file" .fasta)"

        echo "Processing \$file"
        abpoa \$file > "\${base}_consensus.fasta"
        sed -i "1c\\>\${base}" "\${base}_consensus.fasta"
        
    done
    

    for file in *_light_clean.fasta; do
        # Skip if no files matched
        [ -e "\$file" ] || continue

        base="\$(basename "\$file" .fasta)"

        echo "Processing \$file"
        abpoa \$file > "\${base}_consensus.fasta"
        sed -i "1c\\>\${base}" "\${base}_consensus.fasta"
        
    done

    cat *_consensus.fasta > "${meta.well}_all_consensus.fasta"
    """
}