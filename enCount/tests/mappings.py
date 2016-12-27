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
gtfs.get_version_before = Mock(return_value="chM")
gtfs.get_genome_index_dir = Mock(return_value=os.path.join(config.genomes_root, "index", "chM"))
import enCount.externals.junctionseq as junctionseq
junctionseq.sp_call = Mock(return_value=0)


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

        self.sample_name = "SAMP_CHM"
        self.in_fastq_1 = os.path.join(data_root, "fastq", self.sample_name, "ENCFF624OCC.fastq.gz")
        self.in_fastq_2 = os.path.join(data_root, "fastq", self.sample_name, "ENCFF604UQO.fastq.gz")
        self.fastq_pair = (self.in_fastq_1, self.in_fastq_2)
        self.out_mapping_dir = os.path.join(mappings_root, "%s/" % self.sample_name)
        self.out_count_dir = os.path.join(counts_root, "%s/" % self.sample_name)


    def test_count_bam(self):
        # Test counting functionality base on recieved data
        empty_and_create(in_dir=self.out_count_dir)
        in_gtf = gtfs.version_to_path(gtfs.get_version_before(datetime.MINYEAR))
        in_bam = os.path.join(self.out_mapping_dir, mappings.STAR_BAM_NAME)
        mappings.count_bam(in_bam, in_gtf, self.out_count_dir, None)

if __name__ == "__main__":
    unittest.main()
