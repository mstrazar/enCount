from enCount.config import genomes_root
import enCount.gtfs as gtfs
import enCount.externals.dexseq as dexseq

import unittest
import datetime
import os

from mock import Mock
dexseq.sp_call = Mock(return_value=0)
gtfs.get_version_before = Mock(return_value="minimal")

class TestDEXSeq(unittest.TestCase):


    def setUp(self):
        self.gtf_ver = gtfs.get_version_before(datetime.datetime.now())
        self.in_gtf = gtfs.version_to_path(self.gtf_ver)


    def test_gtf_to_gff(self):
        out_gff = os.path.join(genomes_root, "gff", "%s.gff" % self.gtf_ver)
        r = dexseq.gtf_to_gff(self.in_gtf, out_gff)
        self.assertEqual(r, 0)
        self.assertTrue(os.path.exists(out_gff))


