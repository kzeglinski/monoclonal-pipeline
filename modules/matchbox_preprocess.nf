process matchbox_preprocess {
    tag { meta.well }
    label 'process_low'
	publishDir "${params.out_dir}/qc", pattern: "*.tsv",  mode: 'copy', failOnError: true
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://kzeglinski/kzeglinski/nanologix-matchbox:v0.0.4' :
        'ghcr.io/kzeglinski/matchbox:0.1.0' }"

	input:
    tuple val(meta), path(reads)
	val(flanks)
	val(rotate_sequence)
	val(vector_type)
    
    output:
    tuple val(meta), path('*_extracted.fasta'), emit: preprocessed_reads, optional: true
	path("*.tsv"), emit: qc_numbers

	script:
	"""
	# need to edit this to use the flanking sequences from file and rotate sequence
	matchbox \
		--script-file "$projectDir/scripts/${vector_type}_preprocess.mb" \
		-e 0.2 \
		--args "seqid='${meta.well}'" \
		$reads > "${meta.well}_preprocess_qc.tsv" 

	# cat output files (e.g. if heavy and light chain)
	myarray=(`find ./ -maxdepth 1 -name "*.fasta"`)
	if [ \${#myarray[@]} -gt 0 ]; then 
    	cat ./*.fasta > "${meta.well}_extracted.fasta"
	else 
    	echo "no antibody sequences detected"
	fi
	
	"""

}
