import datetime
import sys
import time

import os
import shutil
from mock import Mock

import enCount.config as config
import enCount.db as db
import enCount.encode as encode
import enCount.experiments as experiments
import enCount.gtfs as gtfs
import enCount.integration.experiments_mock as ex_mock
import enCount.mappings as mappings
from enCount.integration.config import sync_test_data

def mock_insert(sample_name="SAMP_CHM", genome_name="chM"):
    """
    Mock insert into db.fastqs and db.gtfs for local testing without download.
    :return:
    """

    # Mock genome index insert
    db.gtfs.insert_one({'gtf_ver': genome_name, 'status': 'ready',
               'time_stamp': datetime.datetime.now(),
               'out_genome_dir': os.path.join(config.genomes_root, "index", genome_name),
               'in_gtf': os.path.join(config.genomes_root, "gtf", "%s.gtf" % genome_name),
               })

    # A hard-coded insert to db.fastqs to prevent download
    db.fastqs.insert_many([
    {'file_path': os.path.join(config.data_root, 'fastq', sample_name, 'ENCFF891EGO.fastq.gz'), 'status': 'ready',
     'size': 2393221685,
     'time_stamp': datetime.datetime(2016, 12, 15, 19, 24, 18, 61080), 'md5': '0ab5a31c32292d4ef03935837401e4d1',
     'url': 'https://www.encodeproject.org/files/ENCFF891EGO/@@download/ENCFF891EGO.fastq.gz', 'e_acc': sample_name,
     'f_acc': 'ENCFF891EGO'},
    {'file_path': os.path.join(config.data_root, 'fastq', sample_name, 'ENCFF667EXM.fastq.gz'), 'status': 'ready',
     'size': 2457367686,
     'time_stamp': datetime.datetime(2016, 12, 15, 19, 24, 18, 83914), 'md5': 'd88f780a73d5d893cc605747593f3d45',
     'url': 'https://www.encodeproject.org/files/ENCFF667EXM/@@download/ENCFF667EXM.fastq.gz', 'e_acc': sample_name,
     'f_acc': 'ENCFF667EXM'},
    {'file_path': os.path.join(config.data_root, 'fastq', sample_name, 'ENCFF513YBI.fastq.gz'), 'status': 'ready',
     'size': 2012751063,
     'time_stamp': datetime.datetime(2016, 12, 15, 19, 24, 18, 85980), 'md5': 'd6e5a7a5295432e95198cf4c474a5122',
     'url': 'https://www.encodeproject.org/files/ENCFF513YBI/@@download/ENCFF513YBI.fastq.gz', 'e_acc': sample_name,
     'f_acc': 'ENCFF513YBI'},
    {'file_path': os.path.join(config.data_root, 'fastq', sample_name, 'ENCFF239VZF.fastq.gz'), 'status': 'ready',
     'size': 2097429588,
     'time_stamp': datetime.datetime(2016, 12, 15, 19, 24, 18, 87278), 'md5': 'fa266ef0de969f7a65371dcecee97718',
     'url': 'https://www.encodeproject.org/files/ENCFF239VZF/@@download/ENCFF239VZF.fastq.gz', 'e_acc': sample_name,
     'f_acc': 'ENCFF239VZF'},
    {'file_path': os.path.join(config.data_root, 'fastq', sample_name, 'ENCFF624OCC.fastq.gz'), 'status': 'ready',
     'size': 1951450497,
     'time_stamp': datetime.datetime(2016, 12, 15, 19, 24, 18, 89384), 'md5': '5d97fc02c972e5fdb9486625d6567523',
     'url': 'https://www.encodeproject.org/files/ENCFF624OCC/@@download/ENCFF624OCC.fastq.gz', 'e_acc': sample_name,
     'f_acc': 'ENCFF624OCC'},
    {'file_path': os.path.join(config.data_root, 'fastq', sample_name, 'ENCFF604UQO.fastq.gz'), 'status': 'ready',
     'size': 2022894136,
     'time_stamp': datetime.datetime(2016, 12, 15, 19, 24, 18, 91356), 'md5': '2cbc7a72d60147a5d3aa196a1b4fdfd6',
     'url': 'https://www.encodeproject.org/files/ENCFF604UQO/@@download/ENCFF604UQO.fastq.gz', 'e_acc': sample_name,
     'f_acc': 'ENCFF604UQO'},
    {'file_path': os.path.join(config.data_root, 'fastq', sample_name, 'ENCFF726LTF.fastq.gz'), 'status': 'ready',
     'size': 2299772019,
     'time_stamp': datetime.datetime(2016, 12, 15, 19, 24, 18, 94333), 'md5': 'bfebdf38bc093f701f09e0ba083a8ef0',
     'url': 'https://www.encodeproject.org/files/ENCFF726LTF/@@download/ENCFF726LTF.fastq.gz', 'e_acc': sample_name,
     'f_acc': 'ENCFF726LTF'},
    {'file_path': os.path.join(config.data_root, 'fastq', sample_name, 'ENCFF569YVH.fastq.gz'), 'status': 'ready',
     'size': 2365716737,
     'time_stamp': datetime.datetime(2016, 12, 15, 19, 24, 18, 96752), 'md5': 'f8b52d8de46ec0c27c9036b1863ce246',
     'url': 'https://www.encodeproject.org/files/ENCFF569YVH/@@download/ENCFF569YVH.fastq.gz', 'e_acc': sample_name,
     'f_acc': 'ENCFF569YVH'}])




class IntExperiments:


    def __init__(self, genome_name="chM", sample_name="SAMP_CHM"):
        """
        Empty databases prior to start
        :return:
        """
        sync_test_data()

        db.fastqs.drop()
        db.mappings.drop()
        db.experiments.drop()
        db.gtfs.drop()

        mock_insert(sample_name, genome_name)
        gtfs.get_version_before = Mock(return_value=genome_name)
        encode.get_online_list = Mock(return_value=ex_mock.get_online_list(sample_name))


        self.e_acc = sample_name
        self.gtf_ver = gtfs.get_version_before(datetime.datetime.now())
        self.design_dir = os.path.join(config.results_root, "junctionseq", self.gtf_ver, self.e_acc)
        self.design_by_sample =  os.path.join(self.design_dir, "decoder.bySample.txt")
        self.design_by_uid = os.path.join(self.design_dir, "decoder.byUID.txt")

        # Design files need to be generated inside the tests
        if os.path.exists(self.design_dir):
            print("Directory %s exists, removing ..." % self.design_dir)
            shutil.rmtree(self.design_dir)


    def test_process(self):
        """
        Fetch list of online experiments, store into database and enqueue for download
        :return:
        """
        # TODO: test when current genome index does not exist

        online_experiments = encode.get_online_list()
        assert len(online_experiments) > 0
        gtf_ver = experiments.add_latest_set(online_experiments)
        assert gtf_ver is not None

        i = 0
        while db.mappings.find({"status": "ready"}).count() < 4:
            print("\n\n\nQueue experiments call %d" % i)
            experiments.process()
            print("\nQueue mappings call %d" % i)
            mappings.process()
            print("Current ready mappings: %d" % db.mappings.find({"status": "ready"}).count())
            time.sleep(2)
            i += 1

        i = 0
        while db.experiments.find({"status": "to process"}).count() > 0:
            print("\n\n\nQueue experiments call %d" % i)
            experiments.process()
            time.sleep(2)
            i += 1

        print("\nRows in database (mappings) after test")
        for row in db.mappings.find({"status": "ready"}):
            print(row)

        print("\nRows in database (experiments) after test")
        for row in db.experiments.find():
            print(row)

        assert os.path.exists(self.design_by_sample)
        assert os.path.exists(self.design_by_uid)


if __name__ == "__main__":
    try:
        genome_name = sys.argv[1]
        sample_name = sys.argv[2]
        test = IntExperiments(sample_name=sample_name, genome_name=genome_name)
    except IndexError:
        print("Genome and sample names not provided, defaulting to chM.")
        test = IntExperiments()

    test.test_process()