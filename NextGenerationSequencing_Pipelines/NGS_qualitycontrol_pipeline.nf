#!/usr/bin/env nextflow

/*
* Author: Sidra S
*/

/*
 * pipeline input parameters
 */
params.reads = "$projectDir/fastq/F3D148_S214_L001_R{1,2}_001.fastq"
params.outdir = "$projectDir/results"
params.fastqcdir = "$projectDir/results/fastqc"
params.cutadaptdir = "$projectDir/results/1_QC-cutadapt"
params.dada2dir = "$projectDir/results/dada2_analysis"



/*
 * Run FastQC
 * 
 */

process FASTQC {
    publishDir params.fastqcdir, mode:'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """
}

process MULTIQC {
    publishDir params.outdir, mode:'copy'

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

process CUTADAPT {
    container 'community.wave.seqera.io/library/pip_cutadapt:66fa80df41437d6b'

    publishDir params.cutadaptdir, mode:'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple path("${sample_id}_trimmed_R1.fastq"), path("${sample_id}_trimmed_R2.fastq")


    script:
    """
    cutadapt -a ^GTGCCAGCMGCCGCGGTAA...ATTAGAWACCCBDGTAGTCC -A ^GGACTACHVGGGTWTCTAAT...TTACCGCGGCKGCTGGCAC -m 100 -M 300 -o  ${sample_id}_trimmed_R1.fastq -p  ${sample_id}_trimmed_R2.fastq ${reads[0]} ${reads[1]}
    """
}

process DADA2 {

    publishDir params.dada2dir, mode:'copy'

    input:
    path(cutadaptdir)
    path(dada2dir)

    output:
    path("ASVs_count.tsv")
    path("summary_tab.txt")
    path("ASVs_AllNichesWithControls.fa")

    script:
    """
    ${projectDir}/amplicon_dada2.R -d ${cutadaptdir} -m ${dada2dir}
    """
}

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    fastqc_ch = FASTQC(read_pairs_ch)
    MULTIQC((fastqc_ch).collect())
    cutadapt_ch = CUTADAPT(read_pairs_ch)

    // read_pairs_ch.view()
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}

