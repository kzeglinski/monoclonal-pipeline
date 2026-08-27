process check_sample_sheet {
    tag "checking sample sheet"
    label 'process_low'

    conda (params.enable_conda ? 'conda-forge::r-tidyverse=2.0.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-tidyverse:2.0.0--33f4d800f6070aac' :
        'community.wave.seqera.io/library/r-tidyverse:2.0.0--dd61b4cbf9e28186' }"

    input:
    path(sample_sheet)

    output:
    path("validated_sample_sheet.csv"), emit: validated_sample_sheet

    script:
    """
    sample_sheet_validation.R ${sample_sheet}
    """
}

process check_flanking_file {
    tag "checking flanking sequences file"
    label 'process_low'

    conda (params.enable_conda ? 'conda-forge::r-tidyverse=2.0.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-tidyverse:2.0.0--33f4d800f6070aac' :
        'community.wave.seqera.io/library/r-tidyverse:2.0.0--dd61b4cbf9e28186' }"

    input:
    path(flanking_sequences)
    val(vector_type)

    output:
    path("flanking_seqs_for_mb.csv"), emit: mb_flanks

    script:
    """
    flanking_file_validation.R ${flanking_sequences} ${vector_type}
    """
}

process concat_reads {
    tag { meta.well }
    label 'process_low'
    publishDir "${params.out_dir}/concat_reads", mode: 'copy', failOnError: true

    input:
    tuple val(meta), path(fastq_files)

    output:
    tuple val(meta), path("${meta.well}.${extn}")

    script:
    if( fastq_files.every { it.name.endsWith('.fastq.gz') } )
        extn = 'fastq.gz'
    else if( fastq_files.every { it.name.endsWith('.fastq') } )
        extn = 'fastq'
    else if( fastq_files.every { it.name.endsWith('.fq.gz') } )
        extn = 'fq.gz'
    else if( fastq_files.every { it.name.endsWith('.fq') } )
        extn = 'fq'
    else
        error "Mixed filetypes or filetypes other than .fastq.gz, .fastq, .fq or .fq.gz is unsupported"
    
    """
    cat ${fastq_files} > "${meta.well}.${extn}"
    """
}

workflow parse_sample_sheet{
// adapted from https://github.com/stevekm/nextflow-demos/blob/master/parse-samplesheet/main.nf
    take:
        sample_sheet
        fastq_dir

    main:
        fastq_extns = [ '.fastq', '.fastq.gz' , '.fq', '.fq.gz' ]
        // validate sample sheet
        check_sample_sheet(sample_sheet)
        validated_sample_sheet = check_sample_sheet.out.validated_sample_sheet

        validated_sample_sheet
            .splitCsv(skip: 1, header: false, sep: ',')
            // or is it better to require exactly matching col names?
            .map{row ->
                def barcode = row[0]
                def well = row[1]
                def sample_id = row[2]
                def vector = row[3]
                // find fastq in the barcoded directories
                def full_path = fastq_dir + "/" + barcode
                def all_files = file(full_path).listFiles()
                def fastq_files = all_files.findAll { fn ->
                        fastq_extns.find { fn.name.endsWith( it ) }
                    }
                return tuple([barcode:barcode, well:well, sample_id:sample_id, vector:vector], fastq_files)
            }
            // concat all files in each barcoded directory
            .set { sample_info } 

        // set vector_type based on whether plate has only antibodies / only nanobodies / mix of both
        // ignore vector types other than antibody and nanobody (i.e. empty_well, irrelevant_phage)
        sample_info
            .map{ it[0].vector }
            .filter{ it == "antibody" || it == "nanobody" }
            .collect()
            .map{ vector_types ->
                if( vector_types.every { it == "antibody" } )
                    return "antibody"
                else if( vector_types.every { it == "nanobody" } )
                    return "nanobody"
                else
                    return "both"
            }
            .set { vector_type }
        
        sample_info_with_reads = concat_reads(sample_info)

    emit:
        sample_info_with_reads
        vector_type
}

workflow parse_flanking_file{
    take:
        flanking_sequences
        vector_type

    main:
        // validate flanking file
        check_flanking_file(flanking_sequences, vector_type)

        flanks = check_flanking_file.out.mb_flanks
            .splitCsv(header: true, sep: ',')
            .map{row -> [(row.parameter): row.sequence]} // convert to map of parameter:sequence
            .reduce{ acc, item -> acc + item } // combine list of maps into one map
        
    emit:
        flanks
}

workflow file_import{
    take:
        sample_sheet
        fastq_dir
        flanking_sequences

    main:
        parsed_sample_sheet = parse_sample_sheet(sample_sheet, fastq_dir)
        samples = parsed_sample_sheet.sample_info_with_reads
        vector_type = parsed_sample_sheet.vector_type
        flanks = parse_flanking_file(flanking_sequences, vector_type)
    
    emit:
        samples // tuple of val(sample metadata), path(concat reads)
        flanks // map of flanking parameter:sequence

}