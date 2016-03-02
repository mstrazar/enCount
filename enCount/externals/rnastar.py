import subprocess
import glob
import os

from .myconfig import NUM_THREADS
from .myconfig import EXECUTABLE_STAR
from .myconfig import GENOME_DIR
from .myconfig import GENOME_FASTA_DIR

def run_star_generate_genome(in_gtf, read_length=100):
    """
    :param in_gtf
        Annotation .gtf file.
    :param read_length
        Read length. Suggested parameter by STAR documentation is 100.

    """
    # List all chromosomes fasta files
    fastas = " ".join(glob.glob(os.path.join(GENOME_FASTA_DIR, "*.fa.gz")))
    assert len(fastas)

    args = [
        "--runThreadN", NUM_THREADS,
        "--runMode", "genomeGenerate",
        "--genomeDir", GENOME_DIR,
        "--genomeFastaFiles", fastas,
        "--sjdbGTFfile", in_gtf,
        "--sjdbOverhang", read_length-1,
    ]

    print(" ".join(args))
    subprocess.call(args)
    return


def run_star(in_fastq, out_dir):
    """
        Run STAR aligner on the in_fastq file.

        Produces out_dir/Aligned.out.sam .

        :param in_fastq
            Input .fastq(.gz) file.
        :param in_gtf

        :param out_dir
            Prefix for the output directory.


    """
    assert in_fastq.endswith(".fastq.gz") or in_fastq.endswith(".fastq")

    args = ["STAR",
            "--readFilesIn", in_fastq,
            "--genomeDir", GENOME_DIR,
            "--runThreadN", NUM_THREADS,
            "--outFileNamePrefix", out_dir,
            "--outSAMtype", "SortedByCoordinate",]

    if in_fastq.endswith(".gz"):
        args.append("--readFilesCommand")
        args.append("zcat")

    print(" ".join(args))
    subprocess.call(args)

    return



