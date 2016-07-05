# coding=utf-8
import subprocess
import enCount.config
import csv
import glob
import os



# FILE STRUCTURE CONSTANTS
EXPERIMENT_ACCESSION = "Experiment accession"
EXPERIMENT_TARGET = "Experiment target"
BIOLOGICAL_REPLICATES = "Biological replicate(s)"
FILE_ACCESSION = "File accession"
CONTROLID = "Non-specific target control-human"


def generate_decoders(in_metafile, control_dir, out_dir):
    """
    :param in_metafile
        Input metadata file containing target (knockdown) experiments
        and controls.
    :param control_dir
        Directory containing metadata files mapping controls to target
        (knockdown) experiments. Assume to be in format
                <ID>.controls.tsv
            e.g.
                ENCSR118EFE.controls.tsv
    :param out_dir
        Output directory

    :result
        Write two JunctionSeq decoder files in out_dir.

        decoder.bySample.txt
            sample.ID       group.ID
            SAMP1   CASE
            SAMP2   CASE
            SAMP3   CTRL
            SAMP4   CTRL
            ...

        sample.ID       lane.ID unique.ID       qc.data.dir     group.ID    input.read.pair.count
            SAMP1   RG1     SAMP1_RG1       SAMP1_RG1       CASE    465298
            SAMP1   RG2     SAMP1_RG2       SAMP1_RG2       CASE    472241
            SAMP2   RG1     SAMP2_RG1       SAMP2_RG1       CASE    461405
            SAMP2   RG2     SAMP2_RG2       SAMP2_RG2       CASE    467713
            ...
    """
    # Unique files containg experimental data in one column
    # Required to merge annotation .gff with nove splices later
    main_unique_ids = set()
    header_bysample = ["sample.ID", "group.ID"]
    header_byuid    = ["sample.ID", "lane.ID", "unique.ID", "qc.data.dir",
                        "group.ID", "input.read.pair.count"]
    out_bysample_main    = open(os.path.join(out_dir, "decoder.bySample.txt"), "wt")
    out_byuid_main       = open(os.path.join(out_dir, "decoder.byUID.txt"), "wt")
    writer_bysample_main = csv.DictWriter(out_bysample_main, delimiter="\t",
                                     fieldnames=header_bysample)
    writer_byuid_main    = csv.DictWriter(out_byuid_main, delimiter="\t",
                                 fieldnames=header_byuid)
    writer_bysample_main.writeheader()
    writer_byuid_main.writeheader()

    # Read unique set of experiment ids
    reader_main = csv.DictReader(open(in_metafile), delimiter="\t")
    experiments = set([row[EXPERIMENT_ACCESSION] for row in reader_main if
                    row[EXPERIMENT_TARGET] != CONTROLID])

    # Individual experimental design files
    for experiment_id in experiments:
        # Search for control ID where experiment_id is contained.
        cid_result = None
        for f in glob.glob(os.path.join(control_dir, "*")):
            cid = os.path.basename(f).split(".")[0]

            reader = csv.DictReader(open(f), delimiter="\t")
            experiments = set([row[EXPERIMENT_ACCESSION] for row in reader])
            if experiment_id in experiments:
                cid_result = cid

        if cid_result is None:
            raise FileNotFoundError("The control does not exist for experiment %s!"
                                    % experiment_id)

        out_dir_exp = os.path.join(out_dir, experiment_id)
        if not os.path.exists(out_dir_exp): os.makedirs(out_dir_exp)
        print("Generating experimental design %s %s" % (experiment_id, out_dir_exp))

        # Create a JunctionSeq decoder bySample file for a given experiment ID
        out_bysample = open(os.path.join(out_dir_exp, "decoder.bySample.txt"), "wt")
        out_byuid = open(os.path.join(out_dir_exp, "decoder.byUID.txt"), "wt")

        writer_bysample = csv.DictWriter(out_bysample, delimiter="\t",
                                         fieldnames=header_bysample)
        writer_byuid = csv.DictWriter(out_byuid, delimiter="\t",
                                     fieldnames=header_byuid)

        writer_bysample.writeheader()
        writer_byuid.writeheader()

        reader = csv.DictReader(open(in_metafile), delimiter="\t")
        for row in sorted(list(reader), key=lambda r: r[EXPERIMENT_TARGET] == CONTROLID):
            if row[EXPERIMENT_ACCESSION] == experiment_id or \
               row[EXPERIMENT_ACCESSION] == cid_result:
                sample_id = row[EXPERIMENT_ACCESSION]
                group_id  = row[EXPERIMENT_TARGET].replace(" ", "_")
                lane_id   = "R"+row[BIOLOGICAL_REPLICATES]
                unique_id = row[FILE_ACCESSION]
                qc_data_dir = row[FILE_ACCESSION]
                input_read_pair_count = "nan"

                row_bysample = { "sample.ID": sample_id, "group.ID": group_id}
                writer_bysample.writerow(row_bysample)

                row_byuid = {"sample.ID": sample_id, "lane.ID":lane_id,
                             "unique.ID": unique_id, "qc.data.dir": qc_data_dir,
                             "group.ID": group_id,
                             "input.read.pair.count": input_read_pair_count}
                writer_byuid.writerow(row_byuid)

                if unique_id not in main_unique_ids:
                    writer_byuid_main.writerow(row_byuid)
                    writer_bysample_main.writerow(row_bysample)
                    main_unique_ids.add(unique_id)

        print("Sucessfully generated %s." % out_bysample)
        print("Sucessfully generated %s." % out_byuid)

    print("Sucessfully generated %s." % out_bysample_main)
    print("Sucessfully generated %s." % out_byuid_main)
    return 0


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
    stdout = open(os.path.join(out_dir, "run_QoRTs_count.out.txt"), "w")
    stderr = open(os.path.join(out_dir, "run_QoRTs_count.err.txt"), "w")

    args = ["java", "-jar", enCount.config.QORTS_JAR,
            "QC", "--stranded",]
    if test_run:
        args = args + ["--testRun"]
    args = args + [in_bam, in_gtf, out_dir]

    print(" ".join(args))
    return subprocess.call(args, stdout=stdout, stderr=stderr)


def run_QoRTs_merge_counts(in_dir, in_decoder, out_dir, test_run=False):
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
    std_dir = os.path.dirname(in_decoder)
    stdout = open(os.path.join(std_dir, "run_QoRTs_merge_counts.out.txt"), "w")
    stderr = open(os.path.join(std_dir, "run_QoRTs_merge_counts.err.txt"), "w")

    args = ["java", "-jar", enCount.config.QORTS_JAR,
            "mergeAllCounts", in_dir, in_decoder, out_dir]
    print(" ".join(args))
    return subprocess.call(args, stdout=stdout, stderr=stderr)


def run_QoRTs_size_factors(in_dir, in_decoder, out_file):
    """
    :param  in_dir
        Directory with raw count files.
    :param in_decoder
        Metadata file indicating sample IDs
    :param out_file
        Size factor file to be generated
    """
    stddir = os.path.dirname(out_file)
    stdout = open(os.path.join(stddir, "run_QoRTs_size_factors.out.txt"), "w")
    stderr = open(os.path.join(stddir, "run_QoRTs_size_factors.err.txt"), "w")

    args = [enCount.config.RSCRIPT, enCount.config.QORTS_R,
            in_dir, in_decoder, out_file]
    print(" ".join(args))
    return subprocess.call(args, stdout=stdout, stderr=stderr)


def run_QoRTs_novel_splices(in_dir, in_gtf, in_size_factors, out_dir,
                            min_count=6):
    """
    Identify novel splice junctions based on minimal read coverage.

    :param in_dir
        Directory with merged count files generated by run_QoRTs_merge_counts.
    :param in_size_factors
        Size factors for all samples. To be used for merging.
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
    stdout = open(os.path.join(out_dir, "run_QoRTs_novel_splices.out.txt"), "w")
    stderr = open(os.path.join(out_dir, "run_QoRTs_novel_splices.err.txt"), "w")

    args = ["java", "-jar", enCount.config.QORTS_JAR, "mergeNovelSplices",
                "--minCount", str(min_count), "--stranded",
                in_dir, in_size_factors, in_gtf, out_dir,]

    print(" ".join(args))
    return subprocess.call(args, stdout=stdout, stderr=stderr)


def run_JunctionSeq_analysis(in_decoder, in_gff, in_gtf_dir, out_dir):
    """
    Run the whole JunctionSeq analysis of differential exon *and* splice
    junction usage.
    :param in_decoder
        Decoder file indicating experimental design.
    :param in_gff.
        .gff file (including novel junctions) as produced by QoRTs.
    :param in_gtf_dir
        Directory with count files with naming as in decoder.
    :param out_dir
        Output directory.
    :result
        Produce .tab files with differential usage analysiss results in the
        output directory.
    """
    stddir = os.path.dirname(in_decoder)
    stdout = open(os.path.join(stddir,"run_JunctionSeq_analysis.out.txt"), "w")
    stderr = open(os.path.join(stddir,"run_JunctionSeq_analysis.err.txt"), "w")

    args = [enCount.config.RSCRIPT, enCount.config.JUNCTIONSEQ_R, in_decoder,
            in_gff, in_gtf_dir, out_dir]

    print(" ".join(args))
    return subprocess.call(args, stdout=stdout, stderr=stderr)

