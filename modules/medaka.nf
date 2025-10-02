// adapted from the nf-core medaka module https://github.com/nf-core/modules/blob/master/modules/nf-core/medaka/main.nf

process medaka {
    tag { meta.well }
    label 'process_medium'
    //publishDir "${params.output_dir}/consensus_sequences", mode: 'copy'


    conda (params.enable_conda ? "bioconda::medaka=1.11.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/abpoa_medaka:dd73dadf53f71452' :
        'community.wave.seqera.io/library/abpoa_medaka:6ebd2179ced25b56' }"

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

        # abpoa consensus starting point
        abpoa \$file > "\${base}_abpoa.fasta"
        sed -i "1c\\>\${base}" "\${base}_abpoa.fasta"

        echo "Processing \$file"
        medaka_consensus -i \$file -d "\${base}_abpoa.fasta" -t $task.cpus -o .
        
        mv consensus.fasta \${base}_consensus.fasta
        
    done
    
    for file in *_light_clean.fasta; do
        # Skip if no files matched
        [ -e "\$file" ] || continue

        base="\$(basename "\$file" .fasta)"

        # abpoa consensus starting point
        abpoa \$file > "\${base}_abpoa.fasta"
        sed -i "1c\\>\${base}" "\${base}_abpoa.fasta"

        echo "Processing \$file"
        medaka_consensus -i \$file -d "\${base}_abpoa.fasta" -t $task.cpus -o .

        mv consensus.fasta \${base}_consensus.fasta
        
    done

    cat *_consensus.fasta > "${meta.well}_all_consensus.fasta"

    """
}