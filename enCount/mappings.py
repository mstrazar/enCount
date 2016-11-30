import os
import enCount
import hashlib
import datetime
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
        print(' problems updating collection mappings record id: {0:s}'.format(dbrec_id))



def map_fastq(fastq_pair, in_genome_dir, out_dir, dbrec_id):

    enCount.externals.rnastar.run_star(
        in_genome_dir=in_genome_dir,
        in_fastq_pair=fastq_pair,
        out_dir=out_dir,
        num_threads=enCount.config.NUM_THREADS)

    new_bam = os.path.join(out_dir, "Aligned.sortedByCoord.out.bam")
    log_final = os.path.join(out_dir, "Log.final.out")

    if os.path.exists(log_final):
        _update_dbrec_status(dbrec_id, 'ready', new_bam)
    else:
        _update_dbrec_status(dbrec_id, 'error', new_bam)
    return

def get_mapping_id(fastq_pair, gtf_ver):
    # create a mapping ID tied to fastq_pair and .gtf version
    base_str = "%s.%s" % (fastq_pair, gtf_ver)
    return hashlib.md5(str(base_str).encode("utf-8")).hexdigest()


def get_bam_file_paths(fastq_pair, gtf_ver):

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
            return mapping['out_dir']
        else:
            # not ready
            return None
    else:
        # Create a new mapping directory
        map_id   = get_mapping_id(fastq_pair, gtf_ver)
        out_dir = os.path.join(enCount.config.mappings_root, gtf_ver, map_id)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        # Check for genome index directory
        gtf_index = enCount.gtfs.get_genome_index_dir(gtf_ver)
        if gtf_index is None:
            # not ready, not inserted
            return None

        time_stamp = datetime.datetime.utcnow()
        new_rec = {'fastq_pair': fastq_pair, 'gtf_ver': gtf_ver, 'gtf_index': gtf_index,
                   'status': 'to map', 'time_stamp': time_stamp,
                   'out_dir': out_dir,
        }
        print('adding new record to mappings collection: {:s}'.format(
            str(new_rec)))
        enCount.db.mappings.insert_one(new_rec)
        # not ready
        return None


def process():
    """Synchronizes database and queue of current mappings jobs."""
    global submitted_mappings
    # query DB to get all records that have status 'to map'
    for e in enCount.db.mappings.find({'status': 'to map'}):
        dbrec_id = str(e['_id'])
        fastq_pair = e['fastq_pair']
        gtf_index = e['gtf_index']
        out_dir = e['out_dir']

        # queue new mappings to process
        if dbrec_id not in submitted_mappings:
            job = enCount.queues.mappings.enqueue_call(
                enCount.mappings.map_fastq,
                args=(fastq_pair, gtf_index, out_dir, dbrec_id),
                result_ttl=-1, ttl=-1, timeout=-1,
            )
            job.meta['mapping_id'] = dbrec_id
            job.save()
            submitted_mappings[dbrec_id] = job

    # clean queue for finished mappings
    for job in enCount.queues.mappings.jobs:
        if job.is_finished or job.is_failed:
            job.cleanup()
            enCount.queues.mappings.remove(job)
