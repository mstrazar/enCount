from enCount.externals.junctionseq import run_QoRTs_count, run_QoRTs_merge, \
    run_QoRTs_size_factors, run_QoRTs_novel_splices

import unittest
from glob import glob
from os import makedirs
from os.path import join, basename, splitext, exists
from shutil import rmtree
from csv import DictReader


class TestQoRTs(unittest.TestCase):
    """ Test QoRTs pipeline """

    data_root = "/Volumes/My Passport/data/encount/"
    data_debug = join(data_root, "externals")
    data_bam = join(data_root, "bam")

    in_decoder_uid = join(data_debug, "QoRTsPipelineWalkthrough", "outputData",
            "JctSeqData", "inst", "extdata", "annoFiles", "decoder.byUID.small.txt")
    in_gtf = join(data_debug,"QoRTsExampleData", "inst", "extdata",
                  "anno.gtf.gz")
    in_dir_bam = data_bam

    data_output = join(data_debug, "output")
    out_dir_raw_cts = join(data_output, "rawCts/")
    out_dir_cts = join(data_output, "cts/")
    out_size_factors = join(data_output, "results", "size_factors.txt")
    out_gtf_dir = join(data_output, "gtf/")

    def setUp(self):
        """ Cleanup test directories """
        for out in [self.out_dir_raw_cts, self.out_dir_cts,
                    self.out_gtf_dir]:
            if exists(out): rmtree(out)
            makedirs(out)

        # Directories must end with /
        assert self.out_dir_raw_cts.endswith("/")
        assert self.out_dir_cts.endswith("/")
        assert self.out_gtf_dir.endswith("/")


    def test_QoRTs(self):
        """ Run counts for single .bam files """
        reader = DictReader(open(self.in_decoder_uid, "rt"),
                            delimiter="\t")
        for row in reader:
            accession = row["unique.ID"]
            in_bam = join(self.in_dir_bam, accession + ".bam")
            out_bam_accession = join(self.out_dir_raw_cts, accession)
            makedirs(out_bam_accession)
            r = run_QoRTs_count(in_bam=in_bam, in_gtf=self.in_gtf,
                            out_dir=out_bam_accession, test_run=True)
            self.assertEqual(r, 0)

        # Merge counts for technical replicates """
        r = run_QoRTs_merge(in_dir=self.out_dir_raw_cts,
                        in_decoder=self.in_decoder_uid,
                        out_dir=self.out_dir_cts)
        self.assertEqual(r, 0)

        # Estimate size factors
        r = run_QoRTs_size_factors(in_dir=self.out_dir_raw_cts,
                              in_decoder=self.in_decoder_uid,
                              out_file=self.out_size_factors)
        self.assertEqual(r, 0)


        # Discover novel splice junctions
        r = run_QoRTs_novel_splices(in_dir=self.out_dir_cts,
                               in_gtf=self.in_gtf,
                               in_size_factors=self.out_size_factors,
                               out_dir=self.out_gtf_dir)
        self.assertEqual(r, 0)

if __name__ == "__main__":
    unittest.main()
