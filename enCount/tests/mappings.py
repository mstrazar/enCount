# coding=utf-8
import os
import enCount.config as config
import enCount.gtfs as gtfs
import enCount.mappings as mappings
import enCount.db as db
from enCount.config import data_root, genomes_root, results_root
import unittest
import datetime
import time

from mock import Mock
gtfs.get_version_before = Mock(return_value="minimal")
gtfs.get_genome_index_dir = Mock(return_value=os.path.join(config.genomes_root, "index", "minimal"))

class TestMappings(unittest.TestCase):
    """
    Test for the gtfs queue and DB.

    Assumes directory structure:
            /endata/genomes
            /endata/genomes/gtf
            /endata/genomes/gtf/
            /endata/genomes/fasta
            /endata/genomes/index
    """

    def setUp(self):
        db.mappings.drop()

        self.in_fastq_1 = os.path.join(data_root, "fastq", "MINIMAL", "ENCFF624OCC_FILTERED.fastq.gz")
        self.in_fastq_2 = os.path.join(data_root, "fastq", "MINIMAL", "ENCFF604UQO_FILTERED.fastq.gz")
        self.fastq_pair = (self.in_fastq_1, self.in_fastq_2)
        self.out_mapping_dir = os.path.join(results_root, "mappings", "MINIMAL/")


    def test_process_queue(self):

        # Insert data on fastqs into database
        in_gtf = gtfs.get_version_before(datetime.MINYEAR)
        out_bam = mappings.get_bam_file_paths(self.fastq_pair, in_gtf)
        self.assertTrue(out_bam is None)

        # Assert correct insertion and directory creation
        lst = list(db.mappings.find({"fastq_pair": self.fastq_pair, "gtf_ver": in_gtf, "status": "to map"}))
        self.assertEqual(len(lst), 1)
        out_bam_dir = lst[0]["out_dir"]
        self.assertTrue(os.path.exists(out_bam_dir))

        # Run mapping jobs
        # TODO: is there a cleaner way to ensure transactions / run in a single container
        mappings.process()
        lst = list(db.mappings.find({"fastq_pair": self.fastq_pair, "gtf_ver": in_gtf, "status": "ready"}))
        while not len(lst):
             time.sleep(1)
             lst = list(db.mappings.find({"fastq_pair": self.fastq_pair, "gtf_ver": in_gtf, "status": "ready"}))

        # Check that results exist
        out_bam_dir = mappings.get_bam_file_paths(self.fastq_pair, in_gtf)
        self.assertTrue(out_bam_dir is not None)




if __name__ == "__main__":
    unittest.main()
