// so that the ascii art etc isn't cluttering up the main script
workflow print_start{
  if(params.help == true){
    log.info """
  ,-^-.   
  |\\/\\|    monoclonal antibody and nanobody sequencing
  `-V-'    analysis pipeline (name tba...)
    H
    H      v0.0.1 (WIP)
    H
 .-;":-.
,'|  `; \

Usage: nextflow run ./main.nf --fastq_dir [input path] --sample_sheet [sample sheet] --antibody_reference [path] --nanobody_reference [path]
--help                 : prints this help message
        
Required arguments:
--out_dir              : where to write output
--fastq_dir            : where fastq files are located (fastq_pass folder)
--sample_sheet         : .csv sample sheet (format: barcode,well,sample_id,vector)
        
Optional
--annotation_method    : use "riot" or "igblast" antibody annotation method (default: riot)
--clone_grouping       : use "vdj_cdr3" or "v_cdr3" to group consensus sequences into clones (default: vdj_cdr3)
--igblast_databases    : location of the igblast databases
--rotate_sequence      : rotate plasmids to begin with this (default: lac promoter)
--flanking_sequences   : .csv flanking file (format: vector_type,flank_name,flank_L,flank_R)
--irrelevant_reference : irrelevant phage reference sequence used for sequencing QC
--antibody_reference   : antibody plasmid reference sequence used for sequencing QC
--nanobody_reference   : nanobody plasmid reference sequence used for sequencing QC
"""
    System.exit(0)
  }

fastq_dir = params.fastq_dir ?: "not provided"
tar_path = params.tar_path ?: "not provided"

log.info """
  ,-^-.   
  |\\/\\|    monoclonal antibody and nanobody sequencing
  `-V-'    analysis pipeline (name tba...)
    H
    H      v0.0.1 (WIP)
    H
 .-;":-.
,'|  `; \

input reads dir     : ${fastq_dir}
input tar file      : ${tar_path}
sample sheet        : ${params.sample_sheet}
antibody reference  : ${params.antibody_reference}
nanobody reference  : ${params.nanobody_reference}
output directory    : ${params.out_dir}
annotation method   : ${params.annotation_method}
clone grouping      : ${params.clone_grouping}
"""

}