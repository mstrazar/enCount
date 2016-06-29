from enCount.externals.junctionseq import *
from enCount.config import data_root

import unittest
from os import makedirs
from os.path import join, exists
from shutil import rmtree
from csv import DictReader


class TestQoRTs(unittest.TestCase):
    """ Test QoRTs pipeline. Testing with provided datasets. """

    # Input
    data_debug = join(data_root, "externals")
    in_dir_bam = join(data_root, "bam")

    # TODO: ensure download of these files environment during Docker build
    in_gtf     = join(data_debug,"QoRTsExampleData", "inst", "extdata",
                  "anno.gtf.gz")
    in_decoder_uid = join(data_debug, "QoRTsPipelineWalkthrough", "outputData",
            "JctSeqData", "inst", "extdata", "annoFiles", "decoder.byUID.small.txt")
    in_decoder_sample = join(data_debug, "QoRTsPipelineWalkthrough",
                             "outputData", "JctSeqData", "inst", "extdata",
                             "annoFiles", "decoder.bySample.small.txt")

    # Output
    data_output = join(data_debug, "output")
    out_dir_raw_cts = join(data_output, "rawCts/")
    out_dir_cts = join(data_output, "cts/")
    out_size_factors = join(data_output, "results", "size_factors.txt")
    out_gtf_dir = join(data_output, "gtf/")
    out_jscs_dir = join(data_output, "jscs/")

    in_count_dir = join(data_output, "gtf/")
    in_gff       = join(data_output, "gtf", "withNovel.forJunctionSeq.gff.gz")


    def setUp(self):
        """ Cleanup test directories """
        for out in [self.out_dir_raw_cts, self.out_dir_cts,
                    self.out_gtf_dir, self.out_jscs_dir]:
            if exists(out): rmtree(out, ignore_errors=True)
            makedirs(out)

        # Directories must end with /
        assert self.in_count_dir.endswith("/")
        assert self.out_dir_raw_cts.endswith("/")
        assert self.out_dir_cts.endswith("/")
        assert self.out_gtf_dir.endswith("/")
        assert self.out_jscs_dir.endswith("/")


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


        # Run JunctionSeq analysis
        r = run_JunctionSeq_analysis(in_count_dir=self.in_count_dir,
                                     in_gff=self.in_gff,
                                     in_decoder=self.in_decoder_sample,
                                     out_dir=self.out_jscs_dir)
        self.assertEqual(r, 0)


if __name__ == "__main__":
    unittest.main()
