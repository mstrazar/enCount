#!/bin/bash

# Test Data
in_gtf="/endata/data/genomes/gtf/Homo_sapiens.GRCh37.75.gtf"
in_fastq_1="/endata/data/fastq/ENCSRSAMPLE/ENCFF018JDZ_filter.fastq.gz"
in_fastq_2="/endata/data/fastq/ENCSRSAMPLE/ENCFF273CYL_filter.fastq.gz"
out_genome_dir="/endata/data/genomes/index/Homo_sapiens.GRCh37.75/"
in_genome_fasta_dir="/endata/data/genomes/fasta/hg19/"
out_dir="/endata/data/bam/"
read_length=100
num_threads=1
clip3pAdapterSeq="-"

# Generate genome index
~/bin/STAR --runThreadN $num_threads --runMode genomeGenerate --genomeDir $out_genome_dir --sjdbGTFfile $in_gtf --sjdbOverhang $(($read_length-1)) --genomeFastaFiles $in_genome_fasta_dir/*.fa

# Run alignment with basic options
~/bin/STAR --readFilesIn $in_fastq_1 $in_fastq_2 --genomeDir $out_genome_dir --runThreadN $num_threads --outFileNamePrefix $out_dir --outSAMtype BAM SortedByCoordinate