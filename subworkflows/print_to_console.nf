// so that the ascii art etc isn't cluttering up the main script
workflow print_start{
  if(params.help == true){
    log.info """
  ,-^-.   
  |\\/\\|    monoclonal nanobody sequencing
  `-V-'    analysis pipeline (name tba...)
    H
    H      v0.0.1 (WIP)
    H
 .-;":-.
,'|  `; \

Usage: nextflow run ./main.nf --fastq_dir [input path] --sample_sheet [sample sheet]
--help                : prints this help message
        
Required arguments:
--out_dir      : where to write output
--fastq_dir    : where fastq files are located (fastq_pass folder)
--sample_sheet : .csv sample sheet (format: barcode,well,sample_id,notes)
--vector_type  : is this an "antibody" or "nanobody" vector
        
Optional (only needed for advanced users)
--igblast_databases   : location of the igblast databases
--rotate_sequence     : rotate plasmids to begin with this (default: lac promoter)
--flanking_sequences  : .csv flanking file (format: vector_type,flank_name,flank_L,flank_R)
--irrelevant_reference: irrelevant reference sequence used for QC
"""
    System.exit(0)
  }

fastq_dir = params.fastq_dir ?: "not provided"
tar_path = params.tar_path ?: "not provided"

log.info """
  ,-^-.   
  |\\/\\|    monoclonal nanobody sequencing
  `-V-'    analysis pipeline (name tba...)
    H
    H      v0.0.1 (WIP)
    H
 .-;":-.
,'|  `; \

input reads dir : ${fastq_dir}
input tar file  : ${tar_path}
sample sheet    : ${params.sample_sheet}
output directory: ${params.out_dir}
vector_type     : ${params.vector_type}
"""

}