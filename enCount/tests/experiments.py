import enCount.encode as encode
import enCount.experiments as experiments
import enCount.tests.experiments_mock as ex_mock
import enCount.db as db

import unittest


from mock import Mock
encode.get_online_list = Mock(return_value=ex_mock.online_experiments)

class TestExperiments(unittest.TestCase):

    def setUp(self):
        # Empty databases prior to start
        db.fastqs.drop()
        db.experiments.drop()

    def test_process(self):
        # Fetch list of online experiments, store into database and enqueue for download
        online_experiments = encode.get_online_list()
        self.assertTrue(len(online_experiments) > 0)
        gtf_ver = experiments.add_latest_set(online_experiments)
        self.assertTrue(gtf_ver is not None)
        experiments.process()