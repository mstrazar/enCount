# coding=utf-8
import os
from enCount.config import data_root, genomes_root, results_root

import unittest
import shutil

from mock import Mock
import enCount.externals.rnastar as rnastar
rnastar.sp_call = Mock(return_value=0)

class TestRNASTAR(unittest.TestCase):

    def setUp(self):
        self.num_threads = 10
        self.in_genome_fasta_dir = os.path.join(genomes_root, "fasta", "minimal")
        self.in_gtf = os.path.join(genomes_root, "gtf", "minimal.gtf")
        self.out_genome_dir = os.path.join(genomes_root, "index", "minimal/")

        self.in_fastq_1 = os.path.join(data_root, "fastq", "MINIMAL", "ENCFF624OCC_FILTERED.fastq.gz")
        self.in_fastq_2 = os.path.join(data_root, "fastq", "MINIMAL", "ENCFF604UQO_FILTERED.fastq.gz")

        self.out_mapping_dir = os.path.join(results_root, "mappings", "MINIMAL/")

        for d in [self.out_genome_dir, self.out_mapping_dir]:
            if not os.path.exists(d):
                os.makedirs(d)

    def test_rnastar_generate_genome(self):

        print("Working directory", os.getcwd())

        # Check input files
        for f_in in [self.in_genome_fasta_dir, self.in_gtf]:
            print("\tChecking %s" % f_in)
            self.assertTrue(os.path.exists(f_in))

        # Generate genome
        r = rnastar.run_star_generate_genome(in_gtf=self.in_gtf,
                                     in_genome_fasta_dir=self.in_genome_fasta_dir,
                                     out_genome_dir=self.out_genome_dir,
                                     num_threads=self.num_threads)
        self.assertEqual(r, 0)



    def test_rnastar_align(self):
        """
        Testing the STAR aligner pipeline inside a container.
        """
        print("Working directory", os.getcwd())

        # Check input files
        for f_in in [self.in_fastq_1, self.in_fastq_2, self.out_genome_dir, self.out_mapping_dir]:
            print("\tChecking %s" % f_in)
            self.assertTrue(os.path.exists(f_in))

        # Run alignment
        r = rnastar.run_star(in_fastq_pair=[self.in_fastq_1, self.in_fastq_2],
                 out_dir=self.out_mapping_dir, in_genome_dir=self.out_genome_dir,
                 num_threads=self.num_threads)

        self.assertEqual(r, 0)

    def tearDown(self):
        for d in [self.out_genome_dir, self.out_mapping_dir]:
            if os.path.exists(d):
                shutil.rmtree(d)

if __name__ == "__main__":
    unittest.main()