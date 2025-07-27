// TO DO: VALIDATE THE SAMPLE SHEET
// TO DO: VALIDATE THE FLANKING SEQUENCE FILE

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
        // validate sample sheet here

        sample_sheet
            .splitCsv(skip: 1, header: false, sep: ',')
            // or is it better to require exactly matching col names?
            .map{row ->
                def barcode = row[0]
                def well = row[1]
                def sample_id = row[2]
                def notes = row[3]
                // find fastq in the barcoded directories
                def full_path = fastq_dir + "/" + barcode
                def all_files = file(full_path).listFiles()
                def fastq_files = all_files.findAll { fn ->
                        fastq_extns.find { fn.name.endsWith( it ) }
                    }
                return tuple([barcode:barcode, well:well, sample_id:sample_id, notes:notes], fastq_files)
            }
            // concat all files in each barcoded directory
            .set { sample_info } 

        sample_info_with_reads = concat_reads(sample_info)

    emit:
        sample_info_with_reads
}

workflow parse_flanking_file{
    take:
        flanking_sequences
        vector_type

    main:
        // validate flanking file here

        flanking_sequences
            .splitCsv(header: true, sep: ',')
            // or is it better to require exactly matching col names?
            .map{row ->
                [
                type: row["vector_type"],
                name: row["flank_name"],
                flank_L: row["flank_L"],
                flank_R: row["flank_R"]
                ]
            }
            .filter{it.type == vector_type}
            .set { flanks }
    emit:
        flanks
}

workflow file_import{
    take:
        sample_sheet
        fastq_dir
        flanking_sequences
        vector_type

    main:
        samples = parse_sample_sheet(sample_sheet, fastq_dir)
        flanks = parse_flanking_file(flanking_sequences, vector_type)
    
    emit:
        samples
        flanks

}