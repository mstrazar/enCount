# coding=utf-8
import sys
import os
import shutil
import unittest
from enCount.config import data_root, genomes_root, results_root

from mock import Mock
import enCount.externals.rnastar as rnastar
rnastar.sp_call = Mock(return_value=0)

class TestRNASTAR(unittest.TestCase):

    genome_lengths = {
        "chM": 16571,
        "ch21ch22": 99434461,
        "initial": 3101804739,
    }

    genome_refs = {
        "chM": 1,
        "ch21ch22": 2,
        "initial": 84,
    }

    counts = {
        "SAMP_CHM": 1000,
        "SAMP_CH21CH22": 5000,
    }

    def setUp(self):
        self.genome_name = "chM"
        self.sample_name = "SAMP_CHM"
        self.in_genome_fasta_dir = os.path.join(genomes_root, "fasta", self.genome_name)
        self.out_genome_dir = os.path.join(genomes_root, "index", self.genome_name)

        self.in_gtf = os.path.join(genomes_root, "gtf", "%s.gtf" % self.genome_name)
        if self.genome_name == "chM": self.in_gtf = None    # contains no junctions, gtf must not be included

        self.in_fastq_1 = os.path.join(data_root, "fastq", self.sample_name, "ENCFF624OCC.fastq.gz")
        self.in_fastq_2 = os.path.join(data_root, "fastq", self.sample_name, "ENCFF604UQO.fastq.gz")

        self.out_mapping_dir = os.path.join(results_root, "mappings", self.sample_name)


    def test_genome_parameters(self):
        """
        Test that correct genome parameters are extracted.
        :return:
        """
        self.assertTrue(os.path.exists(self.in_genome_fasta_dir))
        ln, refs = rnastar._genome_parameters(self.in_genome_fasta_dir)
        self.assertEqual(ln, self.genome_lengths[self.genome_name])
        self.assertEqual(refs, self.genome_refs[self.genome_name])


    def test_rnastar_generate_align(self):
        """
        Testing the STAR aligner pipeline inside a container.
        :return:
        """

        # STEP 1: Generate genome
        print("Working directory", os.getcwd())

        # Check input files
        for f_in in [self.in_genome_fasta_dir]:
            print("\tChecking %s" % f_in)
            self.assertTrue(os.path.exists(f_in))

        # Generate genome ; skip gtf for th chM example with no junctions;
        r = rnastar.run_star_generate_genome(in_gtf=self.in_gtf,
                                             in_genome_fasta_dir=self.in_genome_fasta_dir,
                                             out_genome_dir=self.out_genome_dir)
        self.assertEqual(r, 0)


        # STEP 2: Align reads
        for f_in in [self.in_fastq_1, self.in_fastq_2, self.out_genome_dir, self.out_mapping_dir]:
            print("\tChecking %s" % f_in)
            self.assertTrue(os.path.exists(f_in))

        # Run alignment
        r = rnastar.run_star(in_fastq_pair=[self.in_fastq_1, self.in_fastq_2],
                 out_dir=self.out_mapping_dir, in_genome_dir=self.out_genome_dir)

        self.assertEqual(r, 0)


if __name__ == "__main__":
    unittest.main()