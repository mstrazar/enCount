# coding=utf-8
from enCount.externals.rnastar import run_star, run_star_generate_genome


def test():
    """
    Testing the STAR aligner.
    """

    in_gtf   = "/n/users/martins/Dev/data/encount/QoRTsExampleData/inst/extdata/anno.gtf.gz"
    in_fastq = "/n/users/martins/Dev/data/encount/fastq/ENCFF523XON.fastq"

    out_genome_dir = "/n/users/martins/Dev/data/encode/genomes/"
    in_genome_fasta_dir = "/n/users/martins/Dev/data/genome/hg19_unzipped"


    out_dir   = "/n/users/martins/Dev/data/encount/bam/"

    # Generate genome
    run_star_generate_genome(in_gtf=in_gtf,
                             in_genome_fasta_dir=in_genome_fasta_dir,
                             out_genome_dir=out_genome_dir)

    # Run alignment
    run_star(in_fastq=in_fastq, out_dir=out_dir, in_genome_dir=out_genome_dir)


if __name__ == "__main__":
    test()
