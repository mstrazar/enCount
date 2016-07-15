#!/bin/bash

# Test Data - run from command line
datadir="/endata/"                # in docker
STAR=~/bin/STAR
in_genome_fasta_dir=$datadir/"genomes/fasta/hg19_chrY/"
in_gtf=$datadir/"genomes/gtf/Homo_sapiens.GRCh37.75.chrY.gtf.gz"
in_fastq_1=$datadir/"data/ENCSRSAMPLE/ENCFF018JDZ_filter.fastq.gz"
in_fastq_2=$datadir/"data/ENCSRSAMPLE/ENCFF273CYL_filter.fastq.gz"
out_genome_dir=$datadir/"genomes/index/Homo_sapiens.GRCh37.75/"
read_length=100
num_threads=1
clip3pAdapterSeq="-"


# Generate genome index
$STAR --runThreadN $num_threads --runMode genomeGenerate --genomeDir $out_genome_dir --sjdbGTFfile $in_gtf --sjdbOverhang $(($read_length-1)) --genomeFastaFiles $in_genome_fasta_dir/*.fa

# Run alignment with basic options
out_dir=$datadir/"bam/ENCSRSAMPLE/"
mkdir -p $out_dir
$STAR --readFilesIn $in_fastq_1 $in_fastq_2 --genomeDir $out_genome_dir --runThreadN $num_threads --outFileNamePrefix $out_dir --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat

# Run alignment with complete set of recommended ENCODE options
out_dir=$datadir/"bam/ENCSRSAMPLE-alloptions/"
mkdir -p $out_dir
$STAR --readFilesIn $in_fastq_1 $in_fastq_2 --genomeDir $out_genome_dir --runThreadN $num_threads --outFileNamePrefix $out_dir --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outFilterType BySJout --outFilterMultimapNmax  20 --alignSJoverhangMin  8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000