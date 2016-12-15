# coding=utf-8
import os
import enCount.config as config
import enCount.gtfs as gtfs
import enCount.mappings as mappings
import enCount.db as db
from enCount.config import data_root, mappings_root, counts_root
import unittest
import shutil
import datetime
import time

from mock import Mock
gtfs.get_version_before = Mock(return_value="minimal")
gtfs.get_genome_index_dir = Mock(return_value=os.path.join(config.genomes_root, "index", "minimal"))


def empty_and_create(in_dir):
    # Clear directory
    if os.path.exists(in_dir):
        print("Directory %s exists, emptying ..." % in_dir)
        shutil.rmtree(in_dir)
    os.makedirs(in_dir)


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
        self.in_fastq_1 = os.path.join(data_root, "fastq", "MINIMAL", "ENCFF624OCC.fastq.gz")
        self.in_fastq_2 = os.path.join(data_root, "fastq", "MINIMAL", "ENCFF604UQO.fastq.gz")
        self.fastq_pair = (self.in_fastq_1, self.in_fastq_2)
        self.out_mapping_dir = os.path.join(mappings_root, "MINIMAL/")
        self.out_count_dir = os.path.join(counts_root, "MINIMAL/")



    def test_count_bam(self):
        # Test counting functionality base on recieved data
        empty_and_create(in_dir=self.out_count_dir)
        in_gtf = gtfs.version_to_path(gtfs.get_version_before(datetime.MINYEAR))
        in_bam = os.path.join(self.out_mapping_dir, mappings.STAR_BAM_NAME)
        mappings.count_bam(in_bam, in_gtf, self.out_count_dir, None)
        self.assertTrue(os.path.exists(os.path.join(self.out_count_dir, mappings.QORTS_COUNT_NAME)))

    def test_process_queue(self):

        # Insert data on fastqs into database
        in_gtf = gtfs.get_version_before(datetime.MINYEAR)
        out_bam = mappings.get_bam_file_paths(self.fastq_pair, in_gtf)
        self.assertTrue(out_bam is None)

        # Assert correct insertion and directory creation
        lst = list(db.mappings.find({"fastq_pair": self.fastq_pair, "gtf_ver": in_gtf, "status": "to map"}))
        rec = lst[0]
        self.assertEqual(len(lst), 1)
        out_bam_dir = lst[0]["out_dir"]
        out_count_dir = lst[0]["out_count_dir"]
        empty_and_create(out_bam_dir)
        empty_and_create(out_count_dir)

        # Run mapping and counting jobs by the same process method;
        # first run is map, next run is count
        for upcoming_status in ["to count", "ready"]:
            mappings.process()
            lst = []
            while len(lst) == 0:
                lst = list(db.mappings.find({"fastq_pair": self.fastq_pair,
                                             "gtf_ver": in_gtf, "status": upcoming_status}))
                time.sleep(1)
            print("Current record", lst)

        # Check that results exist and paths are correct
        self.assertTrue(os.path.exists(os.path.join(out_bam_dir, mappings.STAR_BAM_NAME)))
        self.assertTrue(os.path.exists(os.path.join(out_count_dir, mappings.QORTS_COUNT_NAME)))

        # Assert read count is correctly written to database
        rec = lst[0]
        self.assertEqual(rec["read_count"], 1000)


if __name__ == "__main__":
    unittest.main()
