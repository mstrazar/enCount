# coding=utf-8
import os
from enCount.externals.rnastar import run_star, run_star_generate_genome
from enCount.config import data_root, genomes_root, results_root
import shutil
import unittest
import subprocess

class TestRNASTAR(unittest.TestCase):

    def test_rnastar_generate_genome(self):
        in_genome_fasta_dir = os.path.join(genomes_root, "fasta", "hg19")
        in_gtf = os.path.join(genomes_root, "gtf", "initial.gtf")
        out_genome_dir = os.path.join(genomes_root, "index", "initial/")
        num_threads = 10

        print("Working directory", os.getcwd())

        # Check input files
        for f_in in [in_genome_fasta_dir, in_gtf, out_genome_dir]:
            print("\tChecking %s" % f_in)
            self.assertTrue(os.path.exists(f_in))

        # Generate genome
        r = run_star_generate_genome(in_gtf=in_gtf,
                                     in_genome_fasta_dir=in_genome_fasta_dir,
                                     out_genome_dir=out_genome_dir,
                                     num_threads=num_threads)
        self.assertEqual(r, 0)



    def test_rnastar_align(self):
        """
        Testing the STAR aligner pipeline inside a container.
        """
        in_fastq_1 = os.path.join(data_root, "fastq", "ENCSRSAMPLE", "ENCFF018JDZ_filter.fastq.gz")
        in_fastq_2 = os.path.join(data_root, "fastq", "ENCSRSAMPLE", "ENCFF273CYL_filter.fastq.gz")
        in_genome_dir = os.path.join(genomes_root, "index", "initial/")
        out_dir = os.path.join(results_root, "mappings", "ENCSRSAMPLE/")
        num_threads = 10

        print("Working directory", os.getcwd())

        # Check input files
        for f_in in [in_fastq_1, in_fastq_2, in_genome_dir, out_dir]:
            print("\tChecking %s" % f_in)
            self.assertTrue(os.path.exists(f_in))

        # Run alignment
        r = run_star(in_fastq_pair=[in_fastq_1, in_fastq_2],
                 out_dir=out_dir, in_genome_dir=in_genome_dir,
                 num_threads=num_threads)

        self.assertEqual(r, 0)


if __name__ == "__main__":
    unittest.main()