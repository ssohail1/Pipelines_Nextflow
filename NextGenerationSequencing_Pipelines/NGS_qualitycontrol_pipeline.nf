#!/usr/bin/env nextflow

/*
* Author: Sidra S
*/

#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */

// Primary input - text file with all fastq files
params.reads_fastq = "${projectDir}/data/sample_fastq.txt"
params.outdir = "fastqc_results"

/*
 * Generate BAM index file
 */
process run_fastqc {

    container 'community.wave.seqera.io/library/fastqc:0.12.1--af7a5314d5015c29'

    publishDir params.outdir, mode: 'copy'

    input:
        path input_fastq

    output:
        path "${input_fastq}_1_fastqc.html"
        path "${input_fastq}_1_fastqc.zip"
        path "${input_fastq}_2_fastqc.html"
        path "${input_fastq}_2_fastqc.zip"

    script:
    """
    fastqc '$input_fastq' -o '$params.outdir'
    """

}

workflow {

    // Create input channel from a text file listing input file paths
    reads_ch = Channel.fromPath(params.reads_fastq).splitText()

    // Create index file for input BAM file
    run_fastqc(reads_ch)

}
