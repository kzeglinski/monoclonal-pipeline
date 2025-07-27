#!/usr/bin/env nextflow

/*
 * Bring in modules
 */
include { print_start } from './subworkflows/print_to_console'
include { file_import } from './subworkflows/file_import'
include { irrelevant_qc } from './subworkflows/quality_control'
include { seqkit_stats } from './subworkflows/quality_control'
include { matchbox_preprocess } from './modules/matchbox_preprocess'
include { pre_consensus_grouping } from './subworkflows/pre_consensus_grouping'
include { abpoa } from './modules/abpoa'
include { igblast } from './modules/igblast'

/*
 * Run the workflow
 */
workflow {

    // print help message, startup message
    print_start()

    // parameter setup
    // TO-DO: PARAMETER VALIDATION
    out_dir = params.out_dir
    fastq_dir = params.fastq_dir
    sample_sheet = Channel.fromPath(file(params.sample_sheet))
    vector_type = params.vector_type
    igblast_databases = Channel.fromPath(file(params.igblast_databases))

    // collect is used here to turn these into value channels
    // this means we can use the one reference file to process many
    // different fasta inputs
    flanking_sequences = Channel.fromPath(file(params.flanking_sequences)).collect()
    irrelevant_reference = Channel.fromPath(file(params.irrelevant_reference)).collect()
    rotate_sequence = params.rotate_sequence

    // TO-DO: ADD BASECALLING (?)

    // import the files 
    file_import(sample_sheet, fastq_dir, flanking_sequences, vector_type)
    samples = file_import.out.samples
    flanks = file_import.out.flanks

    // quality control
    // basic qc stats 
    all_reads = samples.map{it -> it[1]}.collect()
    qc_stats = seqkit_stats(all_reads)
    irr_phage_qc = irrelevant_qc(samples, irrelevant_reference)

    // read pre-processing
    // TO-DO: EDIT THE PROCESS AND MATCHBOX SCRIPT TO TAKE INTO
    // CONSIDERATION THE FLANKING FILE AND ROTATE SEQUENCE
    // RIGHT NOW IT'S STILL HARD CODED
    preprocessing_results = matchbox_preprocess(samples, flanking_sequences, rotate_sequence, vector_type)
    preprocessing_qc = preprocessing_results.qc_numbers
    preprocessed_reads = preprocessing_results.preprocessed_reads
    
    // qc of the preprocessing and grouping: check if we failed
    // to find ab reads

    // pre-consensus annotation & grouping
    // this is kinda dumb but i can't get it to run an igblast for each
    // well unless i join the channels??
    pp_reads_w_igblast_data = preprocessed_reads.combine(igblast_databases)

    grouped_reads = pre_consensus_grouping(pp_reads_w_igblast_data)
    consensus_input = grouped_reads.consensus_input
    monoclonal_qc = grouped_reads.qc

    // consensus calling
    consensus = abpoa(consensus_input)
    
    // post-consensus annotation
    // just take the combined H/L fasta
    consensus.map{it -> 
        def meta = it[0]
        def comb_cons = it[3]
        return tuple(meta, comb_cons)
    }
    .combine(igblast_databases)
    .set { for_final_annot }

    //for_final_annot.view()

    final_annotation = igblast(for_final_annot, "post")

    // TO-DO: QC REPORT (?)
}

