# coding=utf-8
import os
import enCount.gtfs as gtfs
import enCount.db as db
import enCount.queues as queue
from enCount.config import data_root, genomes_root, results_root
import datetime
import unittest
import shutil
import time


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

    def setUp(self):
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

        # Insert record into database
        self.genome_dir = gtfs.get_genome_index_dir(gtf_ver)
        self.assertTrue(os.path.isdir(self.genome_dir))

        # The method gtfs.get_genome_index_dir is called within gtfs.generate_genome_index
        # and shall not be called elsewhere before a genome index is generated
        gtfs.generate_genome_index(in_gtf, in_genome_fasta_dir, self.genome_dir)
        self.genome_dir = gtfs.get_genome_index_dir(gtf_ver)

        # Check mapping
        mappings = list(db.gtfs.find({"gtf_ver": gtf_ver, "status": "ready"}))
        self.assertEqual(len(mappings), 1)


    def test_process_queue(self):
        """
        Simple test of the genome index generation pipeline using the process queue.
        """
        self.assertTrue(db.gtfs.find().count() == 0)

        # Get latest genome version
        gtf_ver = gtfs.get_version_before(datetime.datetime.min)
        in_gtf = os.path.join(genomes_root, "gtf", "%s.gtf" % gtf_ver)
        self.assertTrue(os.path.exists(in_gtf))

        # Process gtf_version
        self.genome_dir = gtfs.get_genome_index_dir(gtf_ver)
        self.assertTrue(self.genome_dir is not None)
        mappings = list(db.gtfs.find({"gtf_ver": gtf_ver}))
        self.assertEqual(len(mappings), 1)

        # Process outstanding requests; Mock submitted jobs explicitly
        empty = False
        while not empty:
            gtfs.process(mock=True)
            empty = queue.gtfs.is_empty()

        # Wait for database to refresh
        # TODO: is there a cleaner way to ensure transactions?
        mappings = list(db.gtfs.find({"gtf_ver": gtf_ver, "status": "ready"}))
        while not len(mappings):
            time.sleep(1)
            mappings = list(db.gtfs.find({"gtf_ver": gtf_ver, "status": "ready"}))

        # Make sure results exist
        self.genome_dir = gtfs.get_genome_index_dir(gtf_ver)
        self.assertTrue(os.path.isdir(self.genome_dir))


    def tearDown(self):
        shutil.rmtree(self.genome_dir)

if __name__ == "__main__":
    unittest.main()
