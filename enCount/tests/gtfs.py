# coding=utf-8
import os
import enCount.gtfs as gtfs
import enCount.db as db
import enCount.queues as queue
from enCount.config import genomes_root
import datetime
import unittest
import time


# Mock system calls
from mock import Mock
gtfs.rnastar.sp_call = Mock(return_value=0)
gtfs.get_version_before = Mock(return_value="chM")

class TestGtfs(unittest.TestCase):
    """
    Test for the gtfs queue and DB.

    Assumes directory structure:
            /endata/genomes
            /endata/genomes/gtf
            /endata/genomes/gtf/minimal.gtf
            /endata/genomes/fasta
            /endata/genomes/index
    """

    def setUp(self):
        db.gtfs.drop()
        self.gtf_ver = gtfs.get_version_before(datetime.datetime.min)
        self.in_gtf = gtfs.version_to_path(gtf_ver=self.gtf_ver)
        self.in_gff = gtfs.version_to_path(gtf_ver=self.gtf_ver, gff=True)
        self.genome_dir = gtfs.get_genome_index_dir(self.gtf_ver)

        if os.path.exists(self.in_gff):
            os.remove(self.in_gff)


    def test_gtf(self):
        """
        Simple test of the genome index generation pipeline without using the queue mechanism.
        """
        self.assertEqual(db.gtfs.find().count(), 1)
        self.assertTrue(os.path.exists(self.in_gtf))

        # Generate a genome index and get path via DB
        in_genome_fasta_dir = os.path.join(genomes_root, "fasta", self.gtf_ver)
        self.assertTrue(os.path.isdir(in_genome_fasta_dir))

        # Insert record into database
        self.assertTrue(os.path.isdir(self.genome_dir))

        # The method gtfs.get_genome_index_dir is called within gtfs.generate_genome_index
        # and shall not be called elsewhere before a genome index is generated
        gtfs.generate_genome_index(self.in_gtf, in_genome_fasta_dir, self.genome_dir, self.in_gff)
        self.genome_dir = gtfs.get_genome_index_dir(self.gtf_ver)

        # Check mapping
        mappings = list(db.gtfs.find({"gtf_ver": self.gtf_ver, "status": "ready"}))
        self.assertEqual(len(mappings), 1)


    def test_process_queue(self):
        """
        Simple test of the genome index generation pipeline using the process queue.
        """
        self.assertEqual(db.gtfs.find().count(), 1)

        # Get latest genome version
        self.assertTrue(os.path.exists(self.in_gtf))

        # Process gtf_version
        mappings = list(db.gtfs.find({"gtf_ver": self.gtf_ver}))
        self.assertEqual(len(mappings), 1)

        # Process outstanding requests; Mock submitted jobs explicitly
        empty = False
        while not empty:
            gtfs.process(mock=True)
            empty = queue.gtfs.is_empty()

        # Wait for database to refresh
        # TODO: is there a cleaner way to ensure transactions?
        mappings = list(db.gtfs.find({"gtf_ver": self.gtf_ver, "status": "ready"}))
        while not len(mappings):
            time.sleep(1)
            mappings = list(db.gtfs.find({"gtf_ver": self.gtf_ver, "status": "ready"}))

        # Make sure results exist
        self.genome_dir = gtfs.get_genome_index_dir(self.gtf_ver)
        self.assertTrue(os.path.isdir(self.genome_dir))


if __name__ == "__main__":
    unittest.main()
