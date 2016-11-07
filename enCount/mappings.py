import os
import enCount
import hashlib
from bson.objectid import ObjectId

# populate list with currently queued jobs
submitted_mappings = dict(
    (j.meta['mapping_id'], j) for j in enCount.queues.mappings.jobs
)


def _update_dbrec_status(dbrec_id, new_status, new_bam):
    r = enCount.db.mappings.update_one(
        {'_id': ObjectId(dbrec_id)}, {"$set": {'status': new_status,
            'out_bam_file': new_bam}}
    )
    if not r.acknowledged:
        print(' problems updating collection mappings record id: {'
              ':s}'.format(dbrec_id))



def map_fastq(fastq_pair, gtf_ver, out_dir, dbrec_id):

    in_genome_dir = enCount.gtfs.get_genome_index_dir(gtf_ver=gtf_ver)
    enCount.externals.rnastar.run_star(
        in_genome_dir=in_genome_dir,
        in_fastq=fastq_pair,
        out_dir=out_dir,
    )

    new_bam = os.path.join(out_dir, "Aligned.sortedByCoord.out.bam")
    if os.path.exists(new_bam):
        _update_dbrec_status(dbrec_id, 'ready', new_bam)
    else:
        _update_dbrec_status(dbrec_id, 'error', new_bam)
    return


def get_bam_file_path(fastq_pair, gtf_ver):

    """Return path to file or None if file not available."""
    # query DB
    mappings = enCount.db.mappings.find(
        {'fastq_pair': fastq_pair, 'gtf_ver': gtf_ver}
    )
    assert(mappings.count() <= 1)
    mappings = list(mappings)

    if mappings:
        # fetch record from DB
        mapping = mappings[0]
        if mapping['status'] == 'ready':
            return mapping['bam_file_path']
        else:
            # not ready
            return
    else:
        bname   = hashlib.md5(str(fastq_pair)).hexdigest()
        out_dir = os.path.join(enCount.config.mappings_root, gtf_ver, bname)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        time_stamp = datetime.datetime.utcnow()
        new_rec = {'fastq_pair': fastq_pair, 'gtf_ver': gtf_ver,
                   'status': 'to map', 'time_stamp': time_stamp,
                   'out_dir': out_dir,
        }
        print('adding new record to mappings collection: {:s}'.format(
            str(new_rec)))
        enCount.db.mappings.insert_one(new_rec)
        # not ready
        return


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
