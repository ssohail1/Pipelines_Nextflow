#!/usr/bin/env nextflow

/*
 * Bulk RNA-Seq Pipeline
 * Author: Your Name
 * Date: YYYY-MM-DD
 */

nextflow.enable.dsl=2

// Define the workflow
workflow {
    params.reads = file(params.reads)
    params.index = file(params.index)
    params.annotation = file(params.annotation)

    // Step 1: Quality Control
    fastqc(input: params.reads, output: "qc")

    // Step 2: Trimming
    trimmed_reads = trim_adapters(input: params.reads, output: "trimmed")

    // Step 3: Alignment
    aligned_bams = align_reads(input: trimmed_reads, index: params.index, output: "aligned")

    // Step 4: Quantification
    feature_counts(input: aligned_bams, annotation: params.annotation, output: "counts")

    // Step 5: Differential Expression Analysis
    differential_expression(input: "counts", output: "results")
}

// FastQC process
process fastqc {
    publishDir "${output}/qc", mode: 'copy'

    input:
    path reads

    output:
    path "${output}/qc"

    script:
    """
    fastqc -o ${output}/qc ${reads}
    """
}

// Adapter trimming process
process trim_adapters {
    publishDir "${output}/trimmed", mode: 'copy'

    input:
    path reads

    output:
    path "${output}/trimmed"

    script:
    """
    trimmomatic PE -phred33 \
        ${reads[0]} ${reads[1]} \
        ${output}/trimmed/paired_1.fastq.gz ${output}/trimmed/unpaired_1.fastq.gz \
        ${output}/trimmed/paired_2.fastq.gz ${output}/trimmed/unpaired_2.fastq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// Alignment process
process align_reads {
    publishDir "${output}/aligned", mode: 'copy'

    input:
    path reads
    path index

    output:
    path "${output}/aligned/*.bam"

    script:
    """
    hisat2 -x ${index} -1 ${reads[0]} -2 ${reads[1]} -S ${output}/aligned/aligned.sam
    samtools view -bS ${output}/aligned/aligned.sam | samtools sort -o ${output}/aligned/aligned.sorted.bam
    rm ${output}/aligned/aligned.sam
    """
}

// Quantification process
process feature_counts {
    publishDir "${output}/counts", mode: 'copy'

    input:
    path bams
    path annotation

    output:
    path "${output}/counts"

    script:
    """
    featureCounts -a ${annotation} -o ${output}/counts/counts.txt ${bams}
    """
}

// Differential Expression Analysis process (R Script)
process differential_expression {
    publishDir "${output}/results", mode: 'copy'

    input:
    path counts

    output:
    path "${output}/results"

    script:
    """
    Rscript -e '
    library(DESeq2);
    counts <- read.table("${counts}", header=TRUE, row.names=1, comment.char="#");
    coldata <- data.frame(condition=c("control", "control", "treated", "treated"));
    rownames(coldata) <- colnames(counts)[-1];
    dds <- DESeqDataSetFromMatrix(countData=counts[, -1], colData=coldata, design=~condition);
    dds <- DESeq(dds);
    res <- results(dds);
    write.csv(as.data.frame(res), file="${output}/results/differential_expression_results.csv");
    '
    """
}
