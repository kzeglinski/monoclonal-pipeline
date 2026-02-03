#!/usr/bin/env nextflow

/*
 * Bring in modules
 */
include { print_start } from './subworkflows/print_to_console'
include { file_import } from './subworkflows/file_import'
include { irrelevant_qc } from './subworkflows/quality_control'
include { seqkit_stats } from './subworkflows/quality_control'
include { post_consensus_qc } from './subworkflows/quality_control'
include { matchbox_preprocess } from './modules/matchbox_preprocess'
include { preconsensus_group_reads } from './modules/pre_consensus_grouping'
include { abpoa } from './modules/abpoa'
include { igblast } from './modules/igblast'
include { igblast as igblast2 } from './modules/igblast'
include { riot } from './modules/riot'
include { riot as riot2 } from './modules/riot'
include { csv_to_tsv } from './modules/riot'
include { csv_to_tsv as csv_to_tsv2 } from './modules/riot'
include { combine_consensus_seqs } from './modules/cat_outputs'
include { combine_csvs } from './modules/cat_outputs'
include { medaka } from './modules/medaka'

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
    annotation_method = params.annotation_method

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
    if (params.annotation_method == 'riot') {
        initial_annotation_csv = riot(preprocessed_reads, "pre").airr_table
        initial_annotation = csv_to_tsv(initial_annotation_csv).tsv
    }
    else if (params.annotation_method == 'igblast') {
        pp_reads_w_igblast_data = preprocessed_reads.combine(igblast_databases)
        initial_annotation = igblast(pp_reads_w_igblast_data, "pre").airr_table
    }
    else {
        error "Invalid --annotation_method '${params.annotation_method}'. Choose riot or igblast."
    }
    consensus_input = preconsensus_group_reads(initial_annotation).consensus_input

    //monoclonal_qc = grouped_reads.qc
    //monoclonal_qc.map{it -> it[1]}.collect().set{for_combining_csv}
    //combine_csvs(for_combining_csv)

    // consensus calling
    consensus = abpoa(consensus_input)
    //consensus = medaka(consensus_input)

    consensus.map{it -> it[1]}.collect().set{for_combining_fasta}
    combine_consensus_seqs(for_combining_fasta) // combine them all
    
    // post-consensus annotation
    if (params.annotation_method == 'riot') {
        final_annotation_csv = riot2(consensus, "post").airr_table
        final_annotation = csv_to_tsv2(final_annotation_csv).tsv.map{it -> it[1]}.collect()
    }
    else if (params.annotation_method == 'igblast') {
        consensus.combine(igblast_databases).set{ for_final_annot }
        final_annotation = igblast2(for_final_annot, "post").map{it -> it[1]}.collect()
    }
    else {
        error "Invalid --annotation_method '${params.annotation_method}'. Choose riot or igblast."
    }

    // qc of the final annotations
    final_qc = post_consensus_qc(final_annotation)

    // TO-DO: QC REPORT (?)
}

