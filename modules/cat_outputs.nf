// combine outputs as per jing's request

// consensus sequences
process combine_consensus_seqs {
    tag "cat_fastas"
    label 'process_tiny'
    publishDir "${params.out_dir}/consensus", mode: 'copy', pattern: "*_all_sequences.fasta"

    input:
    path(fastas)

    output:
    path("*_all_sequences.fasta")

    script:
    """
    cat *.fasta > 0_all_sequences.fasta
    """
}

// monoclonal 
process combine_csvs {
    tag "cat_csvs"
    label 'process_tiny'
    publishDir "${params.out_dir}/qc", mode: 'copy', pattern: "*_combined_monoclonal_qc.csv"

    input:
    path(csvs)

    output:
    path("*_combined_monoclonal_qc.csv")

    script:
    """
    awk 'FNR==1 && NR!=1{next;}{print}' *.csv > 0_combined_monoclonal_qc.csv
    """
}