# coding=utf-8
import os
from enCount.externals.rnastar import run_star, run_star_generate_genome
from enCount.config import data_root, genomes_root, results_root
import shutil
import unittest

class TestRNASTAR(unittest.TestCase):

    def test_rnastar(self):
        """
        Testing the STAR aligner pipeline inside a container.
        """
        in_genome_fasta_dir = os.path.join(genomes_root, "fasta", "hg19")
        in_gtf = os.path.join(genomes_root, "gtf", "Homo_sapiens.GRCh37.75.gtf")
        in_fastq_1 = os.path.join(data_root, "ENCSRSAMPLE", "ENCFF018JDZ_filter.fastq.gz")
        in_fastq_2 = os.path.join(data_root, "ENCSRSAMPLE", "ENCFF273CYL_filter.fastq.gz")
        out_genome_dir = os.path.join(genomes_root, "index", "Homo_sapiens.GRCh37.75")
        out_dir = os.path.join(results_root, "mappings", "ENCSRSAMPLE")
        num_threads = 10


        print("Working directory", os.getcwd())

        # Check input files
        for f_in in [in_genome_fasta_dir, in_gtf, in_fastq_1, in_fastq_2]:
            self.assertTrue(os.path.exists(f_in))

        # Empty test folders
        for d_out in [out_genome_dir, out_dir]:
            if os.path.exists(d_out):
                shutil.rmtree(d_out)
            os.makedirs(d_out)

        # Generate genome
        r = run_star_generate_genome(in_gtf=in_gtf,
                                 in_genome_fasta_dir=in_genome_fasta_dir,
                                 out_genome_dir=out_genome_dir,
                                 num_threads=num_threads)
        self.assertEqual(r, 0)

        # Run alignment
        r = run_star(in_fastq_pair=[in_fastq_1, in_fastq_2],
                 out_dir=out_dir, in_genome_dir=out_genome_dir,
                 num_threads=num_threads)
        self.assertEqual(r, 0)


if __name__ == "__main__":
    unittest.main()
