# coding=utf-8
import os
from enCount.gtfs import get_version_before, get_genome_index_dir, process
from enCount.config import data_root, genomes_root, results_root
import datetime
import unittest

class TestGtfs(unittest.TestCase):

    def test_gtf(self):
        """
        Testing the genome index generation pipeline.
        """
        # Get initial, minimal version
        gtf_ver = get_version_before(datetime.datetime.min)
        self.assertTrue(os.path.exists(os.path.join(genomes_root, "gtf", "%s.gtf" % gtf_ver)))

        # Create new genome index directory
        get_genome_index_dir(gtf_ver)
        self.assertTrue(os.path.isdir(os.path.join(genomes_root, "index", "initial")))

        # Run process loop
        process()

        # Assert results exist





if __name__ == "__main__":
    unittest.main()
