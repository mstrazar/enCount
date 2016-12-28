# coding=utf-8
import sys
import os
import shutil
from enCount.config import data_root, genomes_root, results_root
import enCount.externals.rnastar as rnastar
from enCount.integration.config import sync_test_data


class IntRNASTAR:

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

    def __init__(self, genome_name="chM", sample_name="SAMP_CHM"):
        sync_test_data()

        self.genome_name = genome_name
        self.sample_name = sample_name
        self.in_genome_fasta_dir = os.path.join(genomes_root, "fasta", self.genome_name)
        self.out_genome_dir = os.path.join(genomes_root, "index", self.genome_name)

        self.in_gtf = os.path.join(genomes_root, "gtf", "%s.gtf" % self.genome_name)
        if self.genome_name == "chM": self.in_gtf = None    # contains no junctions, gtf must not be included

        self.in_fastq_1 = os.path.join(data_root, "fastq", self.sample_name, "ENCFF624OCC.fastq.gz")
        self.in_fastq_2 = os.path.join(data_root, "fastq", self.sample_name, "ENCFF604UQO.fastq.gz")

        self.out_mapping_dir = os.path.join(results_root, "mappings", self.sample_name)

        for d in [self.out_genome_dir, self.out_mapping_dir]:
            if os.path.exists(d):
                print("Removing %s" % d)
                shutil.rmtree(d)
            os.makedirs(d)

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
            assert os.path.exists(f_in)

        # Generate genome ; skip gtf for th chM example with no junctions;
        r = rnastar.run_star_generate_genome(in_gtf=self.in_gtf,
                                             in_genome_fasta_dir=self.in_genome_fasta_dir,
                                             out_genome_dir=self.out_genome_dir)
        assert r == 0


        # STEP 2: Align reads
        for f_in in [self.in_fastq_1, self.in_fastq_2, self.out_genome_dir, self.out_mapping_dir]:
            print("\tChecking %s" % f_in)
            assert os.path.exists(f_in)

        # Run alignment
        r = rnastar.run_star(in_fastq_pair=[self.in_fastq_1, self.in_fastq_2],
                 out_dir=self.out_mapping_dir, in_genome_dir=self.out_genome_dir)
        assert r == 0

        count = rnastar.get_read_count(self.out_mapping_dir)
        print("Number of counted reads: %d" % count)
        assert count == self.counts[self.sample_name]


if __name__ == "__main__":
    try:
        genome_name = sys.argv[1]
        sample_name = sys.argv[2]
        test = IntRNASTAR(sample_name=sample_name, genome_name=genome_name)
    except IndexError:
        print("Genome and sample names not provided, defaulting to chM.")
        test = IntRNASTAR()

    test.test_rnastar_generate_align()