import enCount
import datetime
import os
import csv

# Store control ID for sorting
CONTROLID = "Non-specific target control-human"


def is_updated_set_of_experiments(new_experiments):
    _, cur_experiments = get_latest()

    # check for new experiment ids
    if set(cur_experiments.keys()) != set(new_experiments.keys()):
        return True
    # check for new files or other differences
    for (k, v1) in new_experiments.items():
        v2 = cur_experiments[k]
        if v1 != v2:
            return True
    return False


def get_all():
    """Return list of all experiments and time stamps.

    List is ordered by decreasing time stamp (newest records first).
    """

    # find latest set of experiments (and associated metadata)
    experiments_sets = enCount.db.experiments.find().sort('time_stamp', -1)

    if experiments_sets.count():
        return [(e['time_stamp'], e['experiments']) for e in experiments_sets]
    else:
        return []


def get_latest():
    """Return set of experiments with latest time stamp."""
    all_experiments = get_all()
    if all_experiments:
        return all_experiments[0]
    return None, {}


def add_latest_set(experiments, gtf_ver=None, time_stamp=None):
    """Add new set of experiments obtained from ENCODE.

    Experiments are mapped using the latest gtf available.
    """
    if time_stamp is None:
        time_stamp = datetime.datetime.utcnow()

    if gtf_ver is None:
        gtf_ver = enCount.gtfs.get_version_before(time_stamp)

    enCount.db.experiments.insert_one({
        'time_stamp': time_stamp,
        'experiments': experiments,
        'map_to_gtf': gtf_ver,
        'status': 'to process',
    })

    return gtf_ver


def _get_fastq_files_for_samples(files_recs):
    """

    :param files_recs:
        List of file records retrieved from ENCODE.
        Each record is a line in a metadata file.
    :return:
        A dictionary of biological/technical replicates and associated records.
        Maps samples to fasta files.
    """

    pairings = {}
    for r in files_recs:
        if r['File format'] == 'fastq':
            f_md5 = r['md5sum']
            f_acc = r['File accession']
            f_url = r['File download URL']
            f_size = int(r['Size'])
            b_rep = r['Biological replicate(s)']
            t_rep = r['Technical replicate']
            p_acc = r['Paired with']
            p_end = r['Paired end']
            e_target = r['Experiment target']
            assert p_end == '1' or p_end == '2' or p_end == ''

            rec = (f_acc, f_url, f_size, f_md5, e_target)

            # determine pairing and check consistency
            if p_end == '1':
                k1 = f_acc
                ri = 10         # first mate in pair ; index by file accession
            elif p_end == '2':
                k1 = p_acc
                ri = 11         # second mate in pair ; index by the file in pair
            else:
                k1 = f_acc
                ri = 0          # no pairing ; index by file accession

            # insert into dictionary indexed by
            pairings.setdefault((b_rep, t_rep), {}). \
                     setdefault(k1, []).append((ri, rec))


    # sort paired fastqs
    retd = {}
    for (b_rep, t_rep), tmpd in pairings.items():
        recs = []
        for k1, pairing in tmpd.items():
            pairing.sort()
            assert len(pairing) == 1 or len(pairing) == 2
            if len(pairing) == 1:
                assert pairing[0][0] == 0
            if len(pairing) == 2:
                assert pairing[0][0] == 10 and pairing[1][0] == 11
            recs.append([rec for _, rec in pairing])
        assert retd.setdefault((b_rep, t_rep), recs) == recs
    return retd

# TODO: make sure design files are correctly formatted (uniq. id must be retained). Field names do not match.
def generate_design_files(counts, out_dir):
    """
    Generate an experimental design (decoder) file.

    :param counts
        Input dictionary listing count files and required metadata.
    :param out_dir
        Output directory

    :result
        Write two JunctionSeq decoder files in out_dir.

        sample.ID: Biological replicate
        lane.ID: Technical replicate
        unique.ID: Bio + tech. replicate
        qc.data.dir: Bio + tech. replicate
        group.ID: Experiment target. Make sure control is listed last.
        input.read.pair.count: The number of input reads for this replicate

        decoder.bySample.txt
            sample.ID       group.ID
            SAMP1   CASE
            SAMP1   CASE
            SAMP3   CTRL
            SAMP4   CTRL
            ...

        sample.ID       lane.ID unique.ID       qc.data.dir     group.ID    input.read.pair.count
            SAMP1   RG1     SAMP1_RG1       SAMP1_RG1       CASE    465298
            SAMP1   RG2     SAMP1_RG2       SAMP1_RG2       CASE    472241
            SAMP2   RG1     SAMP2_RG1       SAMP2_RG1       CTRL    461405
            SAMP2   RG2     SAMP2_RG2       SAMP2_RG2       CTRL    467713
            ...
    """
    # Unique files contain experimental data in one column
    out_bysample_main    = open(os.path.join(out_dir, "decoder.bySample.txt"), "wt")
    out_byuid_main       = open(os.path.join(out_dir, "decoder.byUID.txt"), "wt")

    header_bysample = ["sample.ID", "group.ID"]
    header_byuid = ["sample.ID", "lane.ID", "unique.ID", "qc.data.dir", "group.ID", "input.read.pair.count"]
    writer_bysample = csv.DictWriter(out_bysample_main, delimiter="\t", fieldnames=header_bysample)
    writer_byuid    = csv.DictWriter(out_byuid_main, delimiter="\t", fieldnames=header_byuid)
    writer_bysample.writeheader()
    writer_byuid.writeheader()

    sort_key = lambda t: (t[1]["Experiment target"] != CONTROLID,
                         t[1]["Biological replicate"],
                         t[1]["Technical replicate"])

    for cnt_id, cnt_data in sorted(counts.items(), key=sort_key):
        writer_bysample.writerow({"sample.ID": cnt_data["Biological replicate"],
                                  "group.ID": cnt_data["Experiment target"]})

        writer_byuid.writerow({"sample.ID": cnt_data["Biological replicate"],
               "lane.ID": cnt_data["Technical replicate"],
               "unique.ID": "%s_%s" % (cnt_data["Biological replicate"], cnt_data["Technical replicate"]),
               "qc.data.dir": "%s_%s" % (cnt_data["Biological replicate"], cnt_data["Technical replicate"]),
                "group.ID": cnt_data["Experiment target"],
               "input.read.pair.count": cnt_data["Read count"]})

    return 0


def get_design_files(e_acc, e_files, gtf_ver):
    """
    If all BAM and count files are ready,
    a design file can be generated including controls

    :param e_acc:
    :param e_files:
    :param gtf_ver:
    :return:
    """


    fastqs_pairings = _get_fastq_files_for_samples(e_files)
    # From here, all info on e_files is gone (transferred to fastqs_pairings)

    ready = True
    counts = dict()
    for (b_rep, t_rep), pairings in fastqs_pairings.items():
        # Make sure fastq is available, if not will be added
        # to DB and queued for download
        for pairing in pairings:
            fastq_pair = []
            target = None
            for f_acc, f_url, f_size, f_md5, e_target in pairing:
                fastq = enCount.fastqs.get_file_path(e_acc, f_acc, f_url, f_size, f_md5)
                fastq_pair.append(fastq)
                target = e_target

            # All fastqs ready?
            if any((f is None for f in fastq_pair)):
                ready = False
            print("Fastqs ready: %s" % pairings)

            # Get BAM or schedule
            bam = enCount.mappings.get_bam_file_paths(fastq_pair, gtf_ver)
            read_count = enCount.mappings.get_mapping_data(bam, "read_count")
            if bam is None:
                ready = False
            print("BAM ready: %s" % bam)

            # Get count or schedule
            count = enCount.mappings.get_count_file_paths(bam, gtf_ver)
            if count is None:
                ready = False
            print("COUNT ready: %s" % count)

            # All info is available
            counts[count] = {"Biological replicate": b_rep,
                             "Technical replicate": t_rep,
                             "Experiment target": target,
                             "Read count": read_count}

    if not ready:
        return None

    # Create a path to the design file ; rewrite if existing
    # Indexed by gtf experiment accession
    out_dir = os.path.join(enCount.config.results_root, "design", gtf_ver, e_acc)
    if not os.path.exists(out_dir): os.makedirs(out_dir)
    design_file = generate_design_files(counts=counts, out_dir=out_dir)
    return design_file


def process():
    """Synchronizes database and queue of current experiments to process."""
    # query DB to get all records that have status 'to process'
    for e in enCount.db.experiments.find({'status': 'to process'}):

        # Data on 'job' consisting of multiple experiments
        dbrec_id = str(e['_id'])
        experiments = e['experiments']
        gtf_ver = e['map_to_gtf']
        time_stamp = e['time_stamp']

        # collect gtfs that resulted from individual experiments
        gtfs = []
        for e_acc, e_files in experiments.items():

            # Based on the database record, Get a decoder files if all necessary BAM and count files exist

            design_files = get_design_files(e_acc, e_files, gtf_ver)

            # enCount.queues.junctions.enqueue_call()

            # gtf = enCount.junctionseq.get_new_gtf_file_path(bams, gtf_ver)
            # gtfs.append(gtf)

        # merge gtfs if all ready
        if any(gtf is None for gtf in gtfs):
            continue

        # other stuff that needs a new gtf
