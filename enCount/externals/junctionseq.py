# coding=utf-8
import subprocess
import encount.config

def run_QoRTs_count(in_bam, in_gtf, out_dir, test_run=False):
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

    args = ["java", "-jar", encount.config.QORTS_JAR,
            "QC", "--stranded",
            "--testRun" if test_run else "",
            in_bam, in_gtf,
            out_dir]

    print(" ".join(args))
    return subprocess.call(args)


def run_QoRTs_merge(in_dir, in_decoder, out_dir, test_run=False):
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
    args = ["java", "-jar", encount.config.QORTS_JAR,
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
    args = ["/usr/local/bin/Rscript", encount.config.QORTS_R,
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
    args = ["java", "-jar", encount.config.QORTS_JAR, "mergeNovelSplices",
                "--minCount", str(min_count), "--stranded",
                in_dir, in_size_factors, in_gtf, out_dir,]

    print(" ".join(args))
    return subprocess.call(args)


def run_JunctionSeq_analysis(in_decoder, in_gff, in_count_dir, out_dir):
    """
    Run the whole JunctionSeq analysis of differential exon *and* splice
    junction usage.
    :param in_decoder
        Decoder file indicating experimental design.
    :param in_gff.
        .gff file (including novel junctions) as produced by QoRTs.
    :param in_count_dir
        Directory with count files with naming as in decoder.
    :param out_dir
        Output directory.
    :result
        Produce .tab files with differential usage analysiss results in the
        output directory.
    """
    args = ["/usr/local/bin/Rscript", encount.config.JUNCTIONSEQ_R, in_decoder,
            in_gff, in_count_dir, out_dir]

    print(" ".join(args))
    return subprocess.call(args)


