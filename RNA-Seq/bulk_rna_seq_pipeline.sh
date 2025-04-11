#!/bin/bash

# Bulk RNA-Seq Pipeline Script
# Author: Your Name
# Date: YYYY-MM-DD

# Ensure script stops on errors
set -e

# Load necessary modules (if using a cluster with module system)
# module load fastqc hisat2 samtools featurecounts r

# Define directories
RAW_DATA_DIR="raw_data"
TRIMMED_DATA_DIR="trimmed_data"
ALIGNMENT_DIR="alignment"
QC_DIR="quality_control"
COUNTS_DIR="counts"
RESULTS_DIR="results"
REFERENCE_GENOME="reference/genome.fa"
ANNOTATION_FILE="reference/annotation.gtf"

# Create directories
mkdir -p $TRIMMED_DATA_DIR $ALIGNMENT_DIR $QC_DIR $COUNTS_DIR $RESULTS_DIR

# Step 1: Quality Control with FastQC
echo "Running FastQC for raw data..."
fastqc -o $QC_DIR $RAW_DATA_DIR/*.fastq.gz

# Step 2: Trimming with Trimmomatic (optional)
echo "Trimming adapters and low-quality bases..."
for file in $RAW_DATA_DIR/*_R1.fastq.gz; do
  base=$(basename $file _R1.fastq.gz)
  trimmomatic PE -phred33 \
    $RAW_DATA_DIR/${base}_R1.fastq.gz $RAW_DATA_DIR/${base}_R2.fastq.gz \
    $TRIMMED_DATA_DIR/${base}_R1_paired.fastq.gz $TRIMMED_DATA_DIR/${base}_R1_unpaired.fastq.gz \
    $TRIMMED_DATA_DIR/${base}_R2_paired.fastq.gz $TRIMMED_DATA_DIR/${base}_R2_unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

# Step 3: Alignment with HISAT2
echo "Aligning reads to reference genome..."
for file in $TRIMMED_DATA_DIR/*_R1_paired.fastq.gz; do
  base=$(basename $file _R1_paired.fastq.gz)
  hisat2 -x $REFERENCE_GENOME -1 $TRIMMED_DATA_DIR/${base}_R1_paired.fastq.gz -2 $TRIMMED_DATA_DIR/${base}_R2_paired.fastq.gz \
    -S $ALIGNMENT_DIR/${base}.sam
done

# Step 4: Convert SAM to BAM and sort BAM files
echo "Converting SAM to BAM and sorting..."
for file in $ALIGNMENT_DIR/*.sam; do
  base=$(basename $file .sam)
  samtools view -bS $file | samtools sort -o $ALIGNMENT_DIR/${base}.sorted.bam
  rm $file # Remove SAM file to save space
done

# Step 5: Generate read counts with featureCounts
echo "Generating read counts..."
featureCounts -a $ANNOTATION_FILE -o $COUNTS_DIR/counts.txt $ALIGNMENT_DIR/*.sorted.bam

# Step 6: Differential Expression Analysis with R
echo "Running Differential Expression Analysis in R..."
Rscript <<EOF
library(DESeq2)

# Load count data
counts <- read.table("$COUNTS_DIR/counts.txt", header=TRUE, row.names=1, comment.char="#")
coldata <- data.frame(condition=c("control", "control", "treated", "treated"))
rownames(coldata) <- colnames(counts)[-1]

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData=counts[, -1], colData=coldata, design=~condition)

# Run DESeq2 pipeline
dds <- DESeq(dds)
res <- results(dds)

# Save results
write.csv(as.data.frame(res), file="$RESULTS_DIR/differential_expression_results.csv")
EOF

echo "Pipeline completed successfully!"
