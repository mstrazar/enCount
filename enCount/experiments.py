import datetime
import enCount


def is_updated_set_of_experiments(new_experiments):
    cur_rec = get_latest()
    if cur_rec is None:
        cur_experiments = {}
    else:
        _, cur_experiments = latest_rec[0]

    # check if new experiments
    if set(cur_experiments.keys()) != set(new_experiments.keys()):
        return True
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
    """Add new set of experiments obtained from ENCODE."""
    if time_stamp is None:
        time_stamp = datetime.datetime.utcnow()

    enCount.db.experiments.insert_one({
        'time_stamp': time_stamp,
        'experiments': experiments,
    })


def _get_fastq_files(files_recs):
    rets = set()
    for r in files_recs:
        if r['File format'] == 'fastq':
            file_md5 = r['md5sum']
            file_acc = r['File accession']
            file_url = r['File download URL']
            file_size = int(r['Size'])
            rets.add((file_acc, file_url, file_size, file_md5))
    return rets

def process():
    # make sure that all fastq files needed by latest set of experiments are
    # available, if not put them on queue for download
    finished_downloads = set()
    to_remove_because_failed = set()
    for file_name, job in submitted_downloads.items():
        if not job.is_finished:
            continue
        if job.result is not None:
            finished_downloads.add(file_name)
        else:
            # there was some error during downloading
            print('Could not download: {:s}'.format(file_name))
            failed_downloads.add(file_name)
            to_remove_because_failed.add(file_name)

        if enCount.db.fastqs.find({'file_name': file_name},
                                  {'file_name': 1}).limit(1).count() > 0:
            finished_downloads.add(file_name)
    for file_name in finished_downloads:
        job = submitted_downloads.pop(file_name)
        job.cleanup()
        enCount.queues.downloads.remove(job)
        print('Download completed for: {:s}'.format(file_name))
    for file_name in to_remove_because_failed:
        job = submitted_downloads.pop(file_name)
        job.cleanup()
        enCount.queues.downloads.remove(job)

    for e_acc, e_files in latest_experiments.items():
        for f_acc, f_url, f_size, f_md5 in _get_fastq_files(e_files):
            k = (f_acc, f_url, f_size, f_md5)

            target_folder = e_acc
            local_target_folder = os.path.join(enCount.config.data_root,
                                               target_folder)
            # should be relative to enCount.data_root
            if not os.path.isdir(local_target_folder):
                try:
                    os.makedirs(local_target_folder)
                except:
                    print('Error, could not create download folder: '
                          '{:s}'.format(local_target_folder))
                    continue
            target_fname = '{:s}_{:s}_{:d}.fastq.gz'.format(f_acc, f_md5,
                                                            f_size)
            file_name = os.path.join(target_folder, target_fname)

            if file_name in submitted_downloads:
                continue
            if file_name in failed_downloads:
                # do not attempt to download file which failed previously
                continue
            if enCount.db.fastqs.find({'file_name': file_name},
                                      {'file_name': 1}).limit(1).count() == 0:
                print('downloading experiment {:s} file {:s} from '
                      '{:s}'.format(e_acc, f_acc, f_url))
                job = enCount.queues.downloads.enqueue_call(
                    enCount.fastqs.download,
                    args=(f_url, target_folder, target_fname),
                    kwargs={'expected_size': f_size, 'expected_md5': f_md5},
                    result_ttl=-1, ttl=-1, timeout=-1,
                )
                job.meta['file_name'] = file_name
                job.save()
                submitted_downloads[file_name] = job
