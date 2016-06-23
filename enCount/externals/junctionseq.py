# coding=utf-8
import subprocess
from enCount.config import QORTS_JAR, QORTS_R

def run_QoRTs_count(in_bam, in_gtf, out_dir):
    """
    Run QoRTS on aligned .bam files. Functions:
        - Perform QC checks
        - Generate count files

    Require a separate directory for each .bam file acession number.

    :param in_bam
        Input .bam file, one per sample.
    :param in_gtf
        Initial annotation file.
    :param out_dir
        Output directory.

    :result
        Write count.txt file in the output directory.
    """

    args = ["java", "-jar", QORTS_JAR,
            "QC", "--stranded",
            in_bam, in_gtf,
            out_dir]

    print(" ".join(args))
    return subprocess.call(args)


def run_QoRTs_merge(in_dir, in_decoder, out_dir):
    """
    Merge technical replicates from same biological sample.

    :param in_dir
        Directory with raw count files generated by run_QoRTs_count.
    :param in_decoder
        Decoder file experiment name/technical replicate.
    :param out_dir
        Output directory.

    :result
        Merge count files given by the decoder and place in the output
        directory.
    """
    args = ["java", "-jar", QORTS_JAR,
            "mergeAllCounts", in_dir, in_decoder, out_dir]

    print(" ".join(args))
    return subprocess.call(args)


def run_QoRTs_size_factors(in_dir, in_decoder, out_file):
    """
    :param  in_dir
        Directory with raw count files.
    :param in_decoder
        Metadata file indicating sample IDs
    :param out_file
        Size factor file to be generated
    """
    args = ["Rscript", QORTS_R,
            in_dir, in_decoder, out_file]
    print(" ".join(args))
    return subprocess.call(args)


def run_QoRTs_novel_splices(in_dir, in_gtf, in_size_factors, out_dir,
                            min_count=6):
    """
    Identify novel splice junctions based on minimal read coverage.

    :param in_dir
        Directory with merged count files generated by run_QoRTs_merge.
    :param in_size_factors
        Size factors file.
    :param in_gtf
        Initial annotation file.
    :param out_dir
        Output directory.
    :param min_count
        Filtering parameter for determining new junctions, suggested=6.

    :result
        Produce .gff file with novel splice junctions and
        updated count files in the output directory.
    """
    args = ["java", "-jar", QORTS_JAR, "mergeNovelSplices",
                "--minCount", str(min_count), "--stranded",
                in_dir, in_size_factors, in_gtf, out_dir,]

    print(" ".join(args))
    return subprocess.call(args)
