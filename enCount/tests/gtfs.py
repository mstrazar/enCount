# coding=utf-8
import os
import enCount.gtfs as gtfs
import enCount.db as db
from enCount.config import data_root, genomes_root, results_root
import datetime
import unittest
import shutil

# Mock system calls
from mock import Mock
gtfs.rnastar.sp_call = Mock(return_value=0)

class TestGtfs(unittest.TestCase):
    """
    Test for the gtfs queue and DB.

    Assumes directory structure:
            /endata/genomes
            /endata/genomes/gtf
            /endata/genomes/gtf/initial.gtf
            /endata/genomes/fasta
            /endata/genomes/index
    """

    @classmethod
    def setUpClass(cls):
        db.gtfs.drop()

    def test_gtf(self):
        """
        Simple test of the genome index generation pipeline without using the queue mechanism.
        """
        self.assertTrue(db.gtfs.find().count() == 0)

        # Get latest genome version
        gtf_ver = gtfs.get_version_before(datetime.datetime.min)
        in_gtf = os.path.join(genomes_root, "gtf", "%s.gtf" % gtf_ver)
        self.assertTrue(os.path.exists(in_gtf))

        # Generate a genome index and get path via DB
        in_genome_fasta_dir = os.path.join(genomes_root, "fasta")
        self.assertTrue(os.path.isdir(in_genome_fasta_dir))

        # The method gtfs.get_genome_index_dir is called within gtfs.generate_genome_index
        # and shall not be called elsewhere before a genome index is generated
        gtfs.generate_genome_index(in_gtf, in_genome_fasta_dir)
        self.genome_dir = gtfs.get_genome_index_dir(gtf_ver)
        self.assertTrue(os.path.isdir(self.genome_dir))

    def tearDown(self):
        shutil.rmtree(self.genome_dir)

if __name__ == "__main__":
    unittest.main()
