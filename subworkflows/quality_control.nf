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