import os
import time
import signal
import datetime

import enCount

# make sure program gets interrupted in a controlled way
stop_it = False


def signal_handler(signal, frame):
    global stop_it
    print('Ctrl+C pressed. Please wait, will stop...')
    stop_it = True
signal.signal(signal.SIGINT, signal_handler)


# find latest set of experiments (and associated metadata)
experiments_sets = enCount.db.experiments.find().sort('time_stamp', -1)
if experiments_sets.count():
    latest_experiments = experiments_sets[0]['experiments']
    latest_time_stamp = experiments_sets[0]['time_stamp']
    print('Loading latest set of experiments, '
          'with timestamp {:s}'.format(str(latest_time_stamp)))
    print('There are {:d} experiments.'.format(len(latest_experiments)))
else:
    latest_experiments = {}
    print('Warning, no experiment sets present in DB')
    print('Starting from clean state.')
print('')


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


# main scan loop
print('Entering main loop.')
submitted_downloads = dict((j.meta['file_name'], j) for j in
                           enCount.queues.downloads.jobs)
failed_downloads = set()

while not stop_it:
    print('Checking for new data on ENCODE...')
    online_experiments_list = latest_experiments
    #online_experiments_list = enCount.encode.get_online_list()

    # compare latest (in DB) and online (on ENCODE) sets of experiments
    # (including their metadata)
    if enCount.encode.compare_experiments_sets(latest_experiments,
                                               online_experiments_list):
        print('An updated list of experiments is available online.')
        print('Will download all new/changed fastq files needed...')
        print('')
        latest_experiments = online_experiments_list
        enCount.db.experiments.insert_one({
            'time_stamp': datetime.datetime.utcnow(),
            'experiments': latest_experiments,
        })

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

    print('Downloads to process: {:d}'.format(len(submitted_downloads)))
    print('Downloads that failed: {:d}'.format(len(failed_downloads)))
    cn_size = 0
    cn_queued = 0
    cn_started = 0
    cn_finished = 0
    cn_failed = 0
    for j in enCount.queues.downloads.get_jobs():
        cn_size += 1
        cn_queued += j.is_queued
        cn_started += j.is_started
        cn_finished += j.is_finished
        cn_failed += j.is_failed

    print('queue size: {:d}'.format(cn_size))
    print('queued: {:d}'.format(cn_queued))
    print('started: {:d}'.format(cn_started))
    print('finished: {:d}'.format(cn_finished))
    print('failed: {:d}'.format(cn_failed))
    print('')
    time.sleep(10)

print('Stopped.')
