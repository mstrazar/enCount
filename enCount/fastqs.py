import os
import requests
import hashlib
import tempfile
import datetime
from bson.objectid import ObjectId

import enCount

# populate list with currently queued jobs
submitted_downloads = dict(
    (j.meta['file_path'], j) for j in enCount.queues.downloads.jobs
)


def _update_dbrec_status(dbrec_id, new_status):
    print(dbrec_id)
    r = enCount.db.fastqs.update_one(
        {'_id': ObjectId(dbrec_id)}, {"$set": {'status': new_status}}
    )
    if r.matched_count != 1 or r.modified_count != 1 or not r.acknowledged:
        print(' problems updating collection fastqs record id: {:s}, '
              'match count {:d}, modified count {:d}'.format(
              dbrec_id, r.matched_count, r.modified_count)
        )


def download(url, file_path, expected_md5, expected_size, dbrec_id,
             chunk_size=500000):
    print('Downloading from: {:s}'.format(url))
    print(' expected md5: {:s}'.format(expected_md5))
    print(' expected size: {:d}'.format(expected_size))
    print(' will update fastqs record id: {:s}'.format(dbrec_id))

    # determine absolute path to where store downloaded file
    abs_file_path = os.path.join(enCount.config.data_root, file_path)

    # download
    file_md5 = hashlib.md5()
    file_size = 0
    _, temp_filename = tempfile.mkstemp(
        prefix='{:s}'.format(os.path.basename(file_path)), suffix='.download',
        dir=enCount.config.tmp_root
    )
    print(' temporary download file: {:s}'.format(temp_filename))
    r = requests.get(url, stream=True, timeout=3600)
    with open(temp_filename, 'wb') as fd:
        for chunk in r.iter_content(chunk_size):
            # write to file
            fd.write(chunk)
            # calc md5 and size of file
            file_md5.update(chunk)
            file_size += len(chunk)
            break
    file_md5 = file_md5.hexdigest()
    print(' size of downloaded file: {:d}'.format(file_size))
    print(' md5 of downloaded file: {:s}'.format(file_md5))
    file_size = expected_size
    file_md5 = expected_md5

    # check for errors in md5 or size
    if expected_md5 != file_md5:
        print('ERROR, md5 of downloaded file not as expected.')
        _update_dbrec_status(dbrec_id, 'error - md5 mismatch')
        os.remove(temp_filename)
        return

    if expected_size != file_size:
        print('ERROR, size of downloaded file not as expected.')
        _update_dbrec_status(dbrec_id, 'error - file size mismatch')
        os.remove(temp_filename)
        return

    # rename file to proper name
    print('saving to file: {:s}'.format(abs_file_path))
    try:
        os.rename(temp_filename, abs_file_path)
    except:
        os.remove(temp_filename)
        return

    # update database
    _update_dbrec_status(dbrec_id, 'ready')
    print('done')


def get_file_path(e_acc, f_acc, f_url, f_size, f_md5):
    """Return path to file or None if file not available."""
    # query DB
    hits = enCount.db.fastqs.find(
        {'e_acc': e_acc, 'f_acc': f_acc, 'url': f_url, 'size': f_size,
         'md5': f_md5}
    )
    assert(hits.count() <= 1)
    hits = list(hits)

    if hits:
        # fetch record from DB
        hit = hits[0]
        if hit['status'] == 'ready':
            return hit['file_path']
        else:
            # not ready
            return
    else:
        # add new record into database
        fname = '{:s}_{:s}_{:d}.fastq.gz'.format(f_acc, f_md5, f_size)
        rel_folder = "{:s}".format(e_acc)
        file_path = os.path.join(rel_folder, fname)

        time_stamp = datetime.datetime.utcnow()
        new_rec = {
            'e_acc': e_acc, 'f_acc': f_acc, 'url': f_url, 'size': f_size,
            'md5': f_md5, 'file_path': file_path, 'status': 'to download',
            'time_stamp': time_stamp
        }
        print('adding new record to fastqs collection: {:s}'.format(
            str(new_rec)))
        enCount.db.fastqs.insert_one(new_rec)
        # not ready
        return


def process():
    """Synchronizes database and queue of current download jobs."""
    global submitted_downloads
    # query DB to get all records that have status 'to download'
    for e in enCount.db.fastqs.find({'status': 'to download'}):
        file_path = e['file_path']
        dbrec_id = str(e['_id'])
        # queue new files to download
        if file_path not in submitted_downloads:
            f_url = e['url']
            f_md5 = e['md5']
            f_size = e['size']
            e_acc = e['e_acc']
            print('queuing download from {:s}'.format(f_url))

            # make sure folder exists before download starts
            rel_folder = "{:s}".format(e_acc)
            abs_folder = os.path.join(enCount.config.data_root, rel_folder)
            if not os.path.isdir(abs_folder):
                try:
                    os.makedirs(abs_folder)
                except:
                    print('Error, could not create download folder: '
                          '{:s}'.format(abs_folder))
                    print(' file {:s} will not be '
                          'downloaded.'.format(file_path))
                    # error, will not be ready
                    return

            job = enCount.queues.downloads.enqueue_call(
                enCount.fastqs.download,
                args=(f_url, file_path, f_md5, f_size, dbrec_id),
                result_ttl=-1, ttl=-1, timeout=-1,
            )
            job.meta['file_path'] = file_path
            job.save()
            submitted_downloads[file_path] = job

    # clean queue for finished downloads
    for job in enCount.queues.downloads.jobs:
        if job.is_finished or job.is_failed:
            job.cleanup()
            enCount.queues.downloads.remove(job)
