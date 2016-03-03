from enCount.externals.junctionseq import run_QoRTs_count, run_QoRTs_merge, \
    run_QoRTs_size_factors, run_QoRTs_novel_splices
from glob import glob
from os import makedirs
from os.path import join, basename, splitext


def test():

    in_decoder_uid = "/n/users/martins/Dev/data/encount/QoRTsPipelineWalkthrough/outputData/JctSeqData/inst/extdata/annoFiles/decoder.byUID.txt"
    in_dir_bam = "/n/users/martins/Dev/data/encount/bamfiles/"
    in_gtf = "/n/users/martins/Dev/data/encount/QoRTsExampleData/inst/extdata/anno.gtf.gz"

    out_dir_raw_cts = "/n/users/martins/Dev/data/encount/output/rawCts/"
    out_dir_cts = "/n/users/martins/Dev/data/encount/output/cts/"
    out_size_factors = "/n/users/martins/Dev/data/encount/output/results/size_factors.txt"
    out_gtf_dir = "/n/users/martins/Dev/data/encount/output/gtf/"

    # Directories must end with /
    assert  out_dir_raw_cts.endswith("/")
    assert  out_dir_cts.endswith("/")
    assert  out_gtf_dir.endswith("/")


    # Run counts for single .bam files
    for in_bam in glob(join(in_dir_bam, "*.bam")):
        accession         = splitext(basename(in_bam))[0]
        out_bam_accession = join(out_dir_raw_cts, accession)
        makedirs(out_bam_accession)
        run_QoRTs_count(in_bam=in_bam, in_gtf=in_gtf,
                        out_dir=out_bam_accession)

    # Merge counts for technical replicates
    run_QoRTs_merge(in_dir=out_dir_raw_cts,
                    in_decoder=in_decoder_uid,
                    out_dir=out_dir_cts,)

    # Estimate size factors
    run_QoRTs_size_factors(in_dir=out_dir_raw_cts,
                          in_decoder=in_decoder_uid,
                          out_file=out_size_factors)

    # Discover novel splice junctions
    run_QoRTs_novel_splices(in_dir=out_dir_cts,
                           in_gtf=in_gtf,
                           in_size_factors=out_size_factors,
                           out_dir=out_gtf_dir)

if __name__ == "__main__":
    test()
