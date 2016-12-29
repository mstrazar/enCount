# coding=utf-8
import sys
import os
import enCount.config as config
import enCount.gtfs as gtfs
import enCount.mappings as mappings
import enCount.db as db
from enCount.config import data_root, mappings_root, counts_root
from enCount.integration.config import sync_test_data
import shutil
import datetime
import time
from mock import Mock


def empty_and_create(in_dir):
    # Clear directory
    if os.path.exists(in_dir):
        print("Directory %s exists, emptying ..." % in_dir)
        shutil.rmtree(in_dir)
    os.makedirs(in_dir)


class IntMappings:
    """
        Integration tests for mappings.
    """

    def __init__(self, genome_name="chM", sample_name="SAMP_CHM"):
        sync_test_data()
        db.mappings.drop()

        self.sample_name = sample_name
        self.in_fastq_1 = os.path.join(data_root, "fastq", self.sample_name, "ENCFF624OCC.fastq.gz")
        self.in_fastq_2 = os.path.join(data_root, "fastq", self.sample_name, "ENCFF604UQO.fastq.gz")
        self.fastq_pair = (self.in_fastq_1, self.in_fastq_2)
        self.out_mapping_dir = os.path.join(mappings_root, "%s/" % self.sample_name)
        self.out_count_dir = os.path.join(counts_root, "%s/" % self.sample_name)

        gtfs.get_version_before = Mock(return_value=genome_name)
        gtfs.get_genome_index_dir = Mock(return_value=os.path.join(config.genomes_root, "index", genome_name))


    def test_count_bam(self):
        # Test counting functionality base on recieved data
        empty_and_create(in_dir=self.out_count_dir)
        in_gtf = gtfs.version_to_path(gtfs.get_version_before(datetime.MINYEAR))
        in_bam = os.path.join(self.out_mapping_dir, mappings.STAR_BAM_NAME)
        mappings.count_bam(in_bam, in_gtf, self.out_count_dir, None)
        assert os.path.exists(os.path.join(self.out_count_dir, mappings.QORTS_COUNT_NAME))

    def test_process_queue(self):

        # Insert data on fastqs into database
        in_gtf = gtfs.get_version_before(datetime.MINYEAR)
        out_bam = mappings.get_bam_file_paths(self.fastq_pair, in_gtf)
        assert out_bam is None

        # Assert correct insertion and directory creation
        lst = list(db.mappings.find({"fastq_pair": self.fastq_pair, "gtf_ver": in_gtf, "status": "to map"}))
        assert len(lst) == 1
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
        assert os.path.exists(os.path.join(out_bam_dir, mappings.STAR_BAM_NAME))
        assert os.path.exists(os.path.join(out_count_dir, mappings.QORTS_COUNT_NAME))

        # Assert read count is correctly written to database
        rec = lst[0]
        assert rec["read_count"] != -1

if __name__ == "__main__":
    try:
        genome_name = sys.argv[1]
        sample_name = sys.argv[2]
        test = IntMappings(sample_name=sample_name, genome_name=genome_name)
    except IndexError:
        print("Genome and sample names not provided, defaulting to chM.")
        test = IntMappings()

    test.test_count_bam()
    test.test_process_queue()
