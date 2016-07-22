# coding=utf-8
import subprocess
import glob
import os
from enCount.config import STAR_EXEC


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

    tmp_dir = os.path.join(out_genome_dir, "STARtmp")

    args = [STAR_EXEC, "--runThreadN", str(num_threads), "--runMode",
            "genomeGenerate", "--genomeDir", out_genome_dir, "--sjdbGTFfile",
            in_gtf, "--sjdbOverhang", str(read_length-1),
            "--outFileNamePrefix", tmp_dir]

    args.append("--genomeFastaFiles")
    for f in glob.glob(os.path.join(in_genome_fasta_dir, "*.fa")):
        args.append(f)

    print(" ".join(args))
    return subprocess.call(args)


def run_star(in_fastq_pair, in_genome_dir, out_dir, num_threads=4,
             clip3pAdapterSeq="-"):
    """
        Run STAR aligner on the in_fastq file.

        Produces out_dir/Aligned.out.sam .

        :param in_fastq_pair
            Input .fastq file pair.
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
    assert len(in_fastq_pair) == 2
    assert in_fastq_pair[0].endswith(".fastq.gz") or in_fastq_pair[0].endswith(".fastq")
    assert in_fastq_pair[1].endswith(".fastq.gz") or in_fastq_pair[1].endswith(".fastq")

    tmp_dir = os.path.join(out_dir, "STARtmp")

    # Basic options
    args = [STAR_EXEC,
            "--readFilesIn",       in_fastq_pair[0], in_fastq_pair[1],
            "--genomeDir",         in_genome_dir,
            "--runThreadN",        str(num_threads),
            "--outFileNamePrefix", out_dir,
            # "--clip3pAdapterSeq",   clip3pAdapterSeq,
            "--outSAMtype", "BAM", "SortedByCoordinate",
            "--outFileNamePrefix", tmp_dir]

    # Standard ENCODE options (Manual 2.5.1, p. 7)
    args += [
        "--outFilterType", "BySJout",
        "--outFilterMultimapNmax",  "20",
        "--alignSJoverhangMin",  "8",
        "--alignSJDBoverhangMin", "1",
        "--outFilterMismatchNmax", "999",
        "--alignIntronMin", "20",
        "--alignIntronMax", "1000000",
        "--alignMatesGapMax", "1000000",
    ]


    # Process .gzip
    if in_fastq_pair[0].endswith(".gz"):
        args.append("--readFilesCommand")
        args.append("zcat")

    print(" ".join(args))
    return subprocess.call(args)



