#!/usr/bin/env nextflow

/*
 * Bring in modules
 */



/*
 * Run the workflow
 */
workflow {
    /* potentially: basecalling. let's not do it now but ask
     * genomics/whoever if they will do it or if we should roll
     * it in here
     */

    // read in input files and cat them together

    // do some QC to check miniprep/sequencing quality [low priority]

    // with matchbox: rotate, extract VH/VL

    // IgBLAST (pre-consensus)
    
    // Pre-consensus grouping and validation of actual monoclonality

    // consensus calling with abpoa

    // final IgBLAST of consensus sequences (cat them first)

    // report rendering [low priority]
}

