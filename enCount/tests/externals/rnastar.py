# coding=utf-8
import os
import shutil
from enCount.config import data_root, genomes_root, results_root

import unittest
import multiprocessing

from mock import Mock
import enCount.externals.rnastar as rnastar
rnastar.sp_call = Mock(return_value=0)

class TestRNASTAR(unittest.TestCase):

    genome_lengths = {
        "minimal": 16571,
        "initial": 3101804739,
    }

    genome_refs = {
        "minimal": 1,
        "initial": 84,
    }

    def setUp(self):
        self.genome_name = "minimal"
        self.num_threads = multiprocessing.cpu_count()
        self.in_genome_fasta_dir = os.path.join(genomes_root, "fasta", self.genome_name)
        self.out_genome_dir = os.path.join(genomes_root, "index", self.genome_name)

        self.in_gtf = os.path.join(genomes_root, "gtf", "%s.gtf" % self.genome_name)
        if not os.path.exists(self.in_gtf): self.in_gtf = None

        self.in_fastq_1 = os.path.join(data_root, "fastq", "MINIMAL", "ENCFF624OCC.fastq.gz")
        self.in_fastq_2 = os.path.join(data_root, "fastq", "MINIMAL", "ENCFF604UQO.fastq.gz")

        self.out_mapping_dir = os.path.join(results_root, "mappings", "MINIMAL/")

        for d in [self.out_genome_dir, self.out_mapping_dir]:
            if os.path.exists(d) and not isinstance(rnastar.sp_call, Mock):
                print("Removing %s" % d)
                shutil.rmtree(d)
                os.makedirs(d)


    def test_genome_parameters(self):
        self.assertTrue(os.path.exists(self.in_genome_fasta_dir))
        ln, refs = rnastar._genome_parameters(self.in_genome_fasta_dir)
        self.assertEqual(ln, self.genome_lengths[self.genome_name])
        self.assertEqual(refs, self.genome_refs[self.genome_name])



    def test_rnastar_generate_align(self):
        """
        Testing the STAR aligner pipeline inside a container.
        """

        # STEP 1: Generate genome
        print("Working directory", os.getcwd())

        # Check input files
        for f_in in [self.in_genome_fasta_dir]:
            print("\tChecking %s" % f_in)
            self.assertTrue(os.path.exists(f_in))

        # Generate genome ; skip gtf for th minimal example with no junctions;
        r = rnastar.run_star_generate_genome(in_gtf=None,
                                             in_genome_fasta_dir=self.in_genome_fasta_dir,
                                             out_genome_dir=self.out_genome_dir,
                                             num_threads=self.num_threads)
        self.assertEqual(r, 0)


        # STEP 2: Align reads
        for f_in in [self.in_fastq_1, self.in_fastq_2, self.out_genome_dir, self.out_mapping_dir]:
            print("\tChecking %s" % f_in)
            self.assertTrue(os.path.exists(f_in))

        # Run alignment
        r = rnastar.run_star(in_fastq_pair=[self.in_fastq_1, self.in_fastq_2],
                 out_dir=self.out_mapping_dir, in_genome_dir=self.out_genome_dir,
                 num_threads=self.num_threads)

        self.assertEqual(r, 0)

        # Test read count after mapping if not mocked
        if not isinstance(rnastar.sp_call, Mock):
            count = rnastar.get_read_count(self.out_mapping_dir)
            print("Number of counted reads: %d" % count)
            self.assertEqual(count, 1000)

if __name__ == "__main__":
    unittest.main()