import enCount
import datetime
import os
import csv
from bson.objectid import ObjectId

# Store control ID for sorting
CONTROLID = "Non-specific target control-human"
NOVEL_SPLICES = "withNovel.forJunctionSeq.gff.gz"

# Load currently submitted mappings
submitted_mappings = dict(
    (j.meta['mapping_id'], j) for j in enCount.queues.junctions.jobs
)


def _update_dbrec_status(dbrec_id, new_status):
    """
    Update a status in the experiments database.
    :param dbrec_id: Record ID.
    :param new_status: New status: ready, error_sf, error_merge, error_jsed, error.
    :return:
    """
    row = {'status': new_status}

    r = enCount.db.experiments.update_one({'_id': ObjectId(dbrec_id)}, {"$set": row})
    if not r.acknowledged:
        print(' problems updating collection mappings record id: {0:s}'.format(dbrec_id))


def is_updated_set_of_experiments(new_experiments):
    """
    Determine whether new experiments differ from existing experiments in the database.
    :param new_experiments
        Dictionary of new experiments fetched from ENCODE project.
    :return
        Boolean.
    """
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
    """
    Return list of all experiments and time stamps.
    List is ordered by decreasing time stamp (newest records first).
    :return
        List of experiments.
    """

    # find latest set of experiments (and associated metadata)
    experiments_sets = enCount.db.experiments.find().sort('time_stamp', -1)

    if experiments_sets.count():
        return [(e['time_stamp'], e['experiments']) for e in experiments_sets]
    else:
        return []


def get_latest():
    """
    Return set of experiments with latest time stamp.
    :return
        Timestamp, dict of experiments.
    """
    all_experiments = get_all()
    if all_experiments:
        return all_experiments[0]
    return None, {}


def add_latest_set(experiments, gtf_ver=None, time_stamp=None):
    """
    Add new set of experiments obtained from ENCODE.
    Experiments are mapped using the latest gtf available.
    :param experiments
        Dict fo experiments.
    :param gtf_ver
        Gtf version string.
    :param time_stamp
        Time stamp before which the latest gtf is returned.
    :return
        Latest gtf version after insert.
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
    Arrange fastq files into pairs and append other technical information required for processing.
    :param files_recs
        List of file records retrieved from ENCODE.
        Each record is a line in a metadata file.
    :return
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

    Write two JunctionSeq decoder files in out_dir.

        sample.ID: Experiment ID
        lane.ID: Bio. replicate (there are no tech replicates in ENCODE)
        unique.ID: Mapping ID
        qc.data.dir: Mapping (count) dir
        group.ID: Experiment target. Make sure control is listed last.
        input.read.pair.count: The number of input reads for this replicate

        decoder.bySample.txt
            sample.ID       group.ID
            SAMP1   CASE
            SAMP1   CASE
            SAMP2   CTRL
            SAMP2   CTRL
            ...

        sample.ID       lane.ID unique.ID       qc.data.dir     group.ID    input.read.pair.count
            SAMP1   RG1     SAMP1_RG1       SAMP1_RG1       CASE    465298
            SAMP1   RG2     SAMP1_RG2       SAMP1_RG2       CASE    472241
            SAMP2   RG1     SAMP2_RG1       SAMP2_RG1       CTRL    461405
            SAMP2   RG2     SAMP2_RG2       SAMP2_RG2       CTRL    467713
            ...

    :param counts
        Input dictionary listing count files and required metadata.
    :param out_dir
        Output directory
    :return
        Path tuple to design files (byUID, bySample).
    """
    # Unique files contain experimental data in one column

    design_by_sample = os.path.join(out_dir, "decoder.bySample.txt")
    design_by_uid = os.path.join(out_dir, "decoder.byUID.txt")

    out_bysample_main    = open(design_by_sample, "wt")
    out_byuid_main       = open(design_by_uid, "wt")

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

        # Used locally within the experiment folder
        # sample_id = "ctrl" if cnt_data["Experiment target"] == CONTROLID else "case"
        sample_id = cnt_data["File accession"]
        target = cnt_data["Experiment target"].replace(" ", ".")


        writer_bysample.writerow({"sample.ID":  sample_id,
                                  "group.ID":   target,})

        # Tokens should not contain spaces
        # Lane.ID must not be parsed as integer
        writer_byuid.writerow({
               "sample.ID":    sample_id,
               "lane.ID":     "R" + str(cnt_data["Biological replicate"]),
               "unique.ID":   cnt_data["File accession"],
               "qc.data.dir": cnt_data["File accession"],
               "group.ID":    target,
               "input.read.pair.count": cnt_data["Read count"]})


    return design_by_sample, design_by_uid


def get_design_files(e_acc, e_files, gtf_ver):
    """
    If all BAM and count files are ready,
    a design file can be generated including controls

    :param e_acc
        Experiment accession.
    :param e_files
        Experiment files.
    :param gtf_ver
        gtf version string.
    :return
        Path tuple to design files (byUID, bySample).
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
                continue
            print("COUNT ready: %s" % count)

            # All info is available
            counts[count] = {"Experiment accession": e_acc,
                            "File accession": os.path.basename(count),
                            "Biological replicate": b_rep,
                            "Technical replicate": t_rep,
                            "Experiment target": target,
                            "Read count": read_count}

    if not ready:
        return None

    # Create a path to the design file ; rewrite if existing
    # Indexed by gtf experiment accession
    out_dir = os.path.join(enCount.config.results_root, "junctionseq", gtf_ver, e_acc)
    if not os.path.exists(out_dir): os.makedirs(out_dir)
    design_by_sample, design_by_uid = generate_design_files(counts=counts, out_dir=out_dir)
    return design_by_sample, design_by_uid


def QoRTs_pipeline_call(e_acc, gtf_ver, decoder_by_sample, decoder_by_UID, out_jseq_dir, dbrec_id):
    """
    Run QoRTs/JunctionSeq pipeline if the design files are ready. To be enqueued as is.
    Assumes design files are available (not None).
    Updates database record status after finish

    :param e_acc
        Experiment accession.
    :param gtf_ver
        GTF version.
    :param decoder_by_sample
        Design file (by sample).
    :param decoder_by_UID
        Design file (by UID).
    :param out_jseq_dir:
        Output directory to store JunctionSeq results.
    :param dbrec_id
        Database record ID for experiment.
    :return
    """
    # TODO: assert workers see imported modules
    in_count_dir = os.path.join(enCount.config.counts_root, gtf_ver)
    in_gtf = enCount.gtfs.version_to_path(gtf_ver)

    out_dir = os.path.join(enCount.config.results_root, "junctionseq", gtf_ver, e_acc)
    out_size_factors = os.path.join(out_dir, "size_factors.txt")
    out_novel_splices = os.path.join(out_dir, "novel_splices")

    if not os.path.exists(out_novel_splices):
        os.makedirs(out_novel_splices)
    
    # Estimate size factors
    r = enCount.externals.junctionseq.run_QoRTs_size_factors(in_dir=in_count_dir,
                                           in_decoder=decoder_by_UID,
                                           out_file=out_size_factors)
    if r != 0:
        _update_dbrec_status(dbrec_id, "error_sf")
        return

    # Discover novel splice junctions
    s = enCount.externals.junctionseq.run_QoRTs_novel_splices(in_dir=in_count_dir,
                                            in_gtf=in_gtf,
                                            in_size_factors=out_size_factors,
                                            out_dir=out_novel_splices)
    if s != 0:
        _update_dbrec_status(dbrec_id, "error_merge")
        return

    # Run JunctionSeq analysis
    novel_gff = os.path.join(out_novel_splices, NOVEL_SPLICES)
    t = enCount.externals.junctionseq.run_JunctionSeq_analysis(in_count_dir=out_novel_splices,
                                                               in_gff=novel_gff,
                                                               in_decoder=decoder_by_sample,
                                                               out_dir=out_jseq_dir)
    if t != 0:
        _update_dbrec_status(dbrec_id, "error_jseq")
        return

    return


def process():
    """Synchronizes database and queue of current experiments to process."""
    global submitted_mappings

    # query DB to get all records that have status 'to process'
    for e in enCount.db.experiments.find({'status': 'to process'}):

        # Data on 'job' consisting of multiple experiments
        dbrec_id = str(e['_id'])
        experiments = e['experiments']
        gtf_ver = e['map_to_gtf']

        # collect gtfs that resulted from individual experiments
        gtfs = []
        for e_acc, e_files in experiments.items():

            # Based on the database record,
            # get a design files if all necessary BAM and count files exist
            design = get_design_files(e_acc, e_files, gtf_ver)
            if design is None:
                continue

            # TODO: update database entry with design info
            # If design is ready, enqueue a JunctionSeq analysis
            design_by_sample, design_by_uid = design
            design_dir = os.path.dirname(design_by_sample)

            args = (e_acc, gtf_ver, design_by_sample, design_by_uid, design_dir, dbrec_id)

            job_id = dbrec_id
            if job_id not in submitted_mappings:
                job = enCount.queues.junctions.enqueue_call(enCount.experiments.QoRTs_pipeline_call,
                                                      args=args, result_ttl=-1, ttl=-1, timeout=-1)

                job.meta['mapping_id'] = job_id
                job.save()
                submitted_mappings[job_id] = job


        # merge gtfs if all ready
        if any(gtf is None for gtf in gtfs):
            continue

        # other stuff that needs a new gtf
        # ...

    # clean queue for finished mappings
    for job in enCount.queues.junctions.jobs:
        if job.is_finished or job.is_failed:
            job.cleanup()
            enCount.queues.mappings.remove(job)