# coding=utf-8
import subprocess
import glob
import os


def run_star_generate_genome(in_gtf, in_genome_fasta_dir, out_genome_dir,
                             read_length=100, num_threads=4):
    """

    Generate a genome index.

    :param in_gtf
        Annotation .gtf file.
    :param in_genome_fasta_dir
        Directory with genome fasta files.
    :param read_length
        Read length. Suggested parameter by STAR documentation is 100.
    :param num_threads
        Number of threads.
    :param out_genome_dir
        Directory for generating genome indices.

    :results
        Generate genome index files in out_genome_dir.
    """

    # List all chromosomes fasta files
    args = ["STAR",
            "--runThreadN", str(num_threads),
            "--runMode",        "genomeGenerate",
            "--genomeDir",      out_genome_dir,
            "--sjdbGTFfile",    in_gtf,
            "--sjdbOverhang",   str(read_length-1),]

    args.append("--genomeFastaFiles")
    for f in glob.glob(os.path.join(in_genome_fasta_dir, "*.fa")):
        args.append(f)

    print(" ".join(args))
    subprocess.call(args)
    return


def run_star(in_fastq, in_genome_dir, out_dir, num_threads=4,
             clip3pAdapterSeq="-"):
    """
        Run STAR aligner on the in_fastq file.

        Produces out_dir/Aligned.out.sam .

        :param in_fastq
            Input .fastq file.
        :param in_genome_dir
            Directory with generated genome indices.
        :param num_threads
            Number of threads.
        :param out_dir
            Prefix for the output directory.

        :param clip3pAdapterSeq
            string(s): adapter sequences to clip from 3p of each mate.
            If one value is given, it will be assumed the same for both mates.
            Default: -

        :result
            Generate a .bam file sorted by coordinate in out_dir.
            Assume 3' adaptor clipping.

    """
    assert in_fastq.endswith(".fastq.gz") or in_fastq.endswith(".fastq")

    args = ["STAR",
            "--readFilesIn",       in_fastq,
            "--genomeDir",         in_genome_dir,
            "--runThreadN",        str(num_threads),
            "--outFileNamePrefix", out_dir,
            "--clip3pAdapterSeq",   clip3pAdapterSeq,
            "--outSAMtype", "BAM", "SortedByCoordinate",]

    if in_fastq.endswith(".gz"):
        args.append("--readFilesCommand")
        args.append("zcat")

    print(" ".join(args))
    subprocess.call(args)

    return



