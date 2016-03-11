import os

import enCount

# populate list with currently queued jobs
submitted_mappings = dict(
    (j.meta['mapping_id'], j) for j in enCount.queues.mappings.jobs
)


def _update_dbrec_status(dbrec_id, new_status):
    r = enCount.db.mappings.update_one(
        {'_id': ObjectId(dbrec_id)}, {"$set": {'status': new_status}}
    )
    if not r.acknowledged:
        print(' problems updating collection mappings record id: {'
              ':s}'.format(dbrec_id))

def map_fastq(fastqs, gtf_ver, outbam_fname):
    return


def get_bam_file_path(fastqs, gtf_ver):
    return


def process():
    """Synchronizes database and queue of current mappings jobs."""
    global submitted_mappings
    # query DB to get all records that have status 'to map'
    for e in enCount.db.mappings.find({'status': 'to map'}):
        mapping_id = str(e['_id'])
        fastqs = e['fastqs']
        gtf_ver = e['gtf_ver']
        outbam_fname = e['outbam_fname']

        # queue new mappings to process
        if mapping_id not in submitted_mappings:
            e_acc = e['e_acc']

            # make sure folder exists before download starts
            rel_folder = os.path.join("mappings", gtf_ver)
            abs_folder = os.path.join(enCount.config.results_root, rel_folder)
            outbam_fname = \
                os.path.join(rel_folder, '{:s}.bam'.format(mapping_id))

            if not os.path.isdir(abs_folder):
                try:
                    os.makedirs(abs_folder)
                except:
                    print('Error, could not create mapping folder: '
                          '{:s}'.format(abs_folder))
                    print(' file {:s} will not be '
                          'generated.'.format(outbam_fname))
                    # error, will not be ready
                    return

            job = enCount.queues.enqueue_call(
                enCount.mappings.map_fastq,
                args=(fastqs, gtf_ver, outbam_fname),
                result_ttl=-1, ttl=-1, timeout=-1,
            )
            job.meta['mapping_id'] = mapping_id
            job.save()
            submitted_mappings[mapping_id] = job

    # clean queue for finished mappings
    for job in enCount.queues.mappings.jobs:
        if job.is_finished or job.is_failed:
            job.cleanup()
            enCount.queues.mappings.remove(job)
