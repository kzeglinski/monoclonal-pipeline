process matchbox_preprocess {
    tag { meta.well }
    label 'process_low'
	publishDir "${params.out_dir}/sequencing_qc", pattern: "*.tsv",  mode: 'copy', failOnError: true
    container "ghcr.io/kzeglinski/matchbox:0.2.0"

	input:
    tuple val(meta), path(reads)
	val(flanks)
	// val(rotate_sequence)
	// val(vector_type)
    
    output:
    tuple val(meta), path('*_extracted.fasta'), emit: preprocessed_reads, optional: true
	path("*.tsv"), emit: qc_numbers

	script:
	"""
	# add rotate sequence?

	if [ "${meta.vector}" == "nanobody" ]; then
		matchbox \
			--script-file "$projectDir/scripts/nanobody_preprocess.mb" \
			-e 0.25 \
			--args "seqid='${meta.well}', nb_flank_L='${flanks.nanobody_nb_L}', nb_flank_R='${flanks.nanobody_nb_R}'" \
			$reads > "${meta.well}_preprocess_qc.tsv"
	elif [ "${meta.vector}" == "antibody" ]; then
		matchbox \
			--script-file "$projectDir/scripts/antibody_preprocess.mb" \
			-e 0.25 \
			--args "seqid='${meta.well}', vh_flank_L='${flanks.antibody_vh_L}', vh_flank_R='${flanks.antibody_vh_R}', vlk_flank_L='${flanks.antibody_vlk_L}', vlk_flank_R='${flanks.antibody_vlk_R}', vll_flank_L='${flanks.antibody_vll_L}', vll_flank_R='${flanks.antibody_vll_R}'" \
			$reads > "${meta.well}_preprocess_qc.tsv"
	else
		echo "vector type is ${meta.vector}; skipping matchbox preprocessing" > "${meta.well}_preprocess_qc.tsv"
	fi

	# cat output fasta files
	myarray=(`find ./ -maxdepth 1 -name "*.fasta"`)
	if [ \${#myarray[@]} -gt 0 ]; then 
    	cat ./*.fasta > "${meta.well}_extracted.fasta"
	fi
	
	"""

}
