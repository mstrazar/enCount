import datetime
import enCount


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


def add_latest_set(experiments, time_stamp=None):
    """Add new set of experiments obtained from ENCODE.

    Experiments are mapped using the latest gtf available.
    """
    if time_stamp is None:
        time_stamp = datetime.datetime.utcnow()

    gtf_ver = enCount.gtfs.get_version_before(time_stamp)
    enCount.db.experiments.insert_one({
        'time_stamp': time_stamp,
        'experiments': experiments,
        'map_to_gtf': gtf_ver,
        'status': 'to process',
    })


def _get_fastq_files_for_samples(files_recs):
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
            assert p_end == '1' or p_end == '2' or p_end == ''

            rec = (f_acc, f_url, f_size, f_md5)

            # determine pairing and check consistency
            if p_end == '1':
                k1 = f_acc
                ri = 10
            elif p_end == '2':
                k1 = p_acc
                ri = 11
            else:
                k1 = f_acc
                ri = 0
            pairings.setdefault((b_rep, t_rep), {}). \
                     setdefault(k1, []).append((ri, rec))

    # sort paired fastqs
    retd = {}
    for (s, t_rep), tmpd in pairings.items():
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


def process():
    """Synchronizes database and queue of current experiments to process."""
    # query DB to get all records that have status 'to process'
    for e in enCount.db.experiments.find({'status': 'to process'}):
        dbrec_id = str(e['_id'])
        experiments = e['experiments']
        gtf_ver = e['map_to_gtf']
        time_stamp = e['time_stamp']

        # collect gtfs that resulted from individual experiments
        gtfs = []
        for e_acc, e_files in experiments.items():
            bams = []
            fastqs_pairings = _get_fastq_files_for_samples(e_files)
            for (b_rep, t_rep), pairings in fastqs_pairings.items():
                # make sure fastq is available, if not will be added
                # to DB and queued for download
                for pairing in pairings:
                    fastqs = []
                    for f_acc, f_url, f_size, f_md5 in pairing:
                        fastq = enCount.fastqs.get_file_path(e_acc, f_acc,
                                                             f_url, f_size,
                                                             f_md5)
                        fastqs.append(fastq)

                    # map only if all fastqs ready
                    if any((f is None for f in fastqs)):
                        continue
                    bam = enCount.mappings.get_bam_file_path(fastqs, gtf_ver)
                    bams.append((b_rep, t_rep, bam))

            # run junction if all bam files ready
            if any(bam is None for _, _, bam in bams):
                continue
            # gtf = enCount.junctionseq.get_new_gtf_file_path(bams, gtf_ver)
            # gtfs.append(gtf)

        # merge gtfs if all ready
        if any(gtf is None for gtf in gtfs):
            continue

        # other stuff that needs a new gtf
