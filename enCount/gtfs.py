import os
import enCount
import datetime
import time
from bson.objectid import ObjectId
from enCount.config import genomes_root
import sys

# Externals
from enCount.externals import rnastar

submitted_gtf_generates = dict(
    (j.meta['genome_index_id'], j) for j in enCount.queues.gtfs.jobs
)

def _update_dbrec_status(dbrec_id, new_status):
    """
    Update a database record on a .gtf file.
    """
    r = enCount.db.gtfs.update_one(
        {'_id': ObjectId(dbrec_id)}, {"$set": {'status': new_status}}
    )
    if r.acknowledged:
        print("Updated object in collection gtfs: %s with status %s" % (str(dbrec_id), new_status))
    else:
        print('Problems updating collection gtfs record id: {0:s}'.format(dbrec_id))


def generate_genome_index(in_gtf, in_genome_fasta_dir, out_genome_dir, mock=False):
    """
    Generate a genome index using STAR for a given .gtf file.
    """
    # Find DB record ID which must exist prior to method call
    mappings = list(enCount.db.gtfs.find({"in_gtf": in_gtf}))
    assert len(mappings) == 1
    dbrec_id = mappings[0]["_id"]

    if mock:
        r = 0
    else:
        r = rnastar.run_star_generate_genome(in_gtf=in_gtf,
                                 in_genome_fasta_dir=in_genome_fasta_dir,
                                 out_genome_dir=out_genome_dir,
                                 num_threads=enCount.config.NUM_THREADS)

    if r == 0:
        _update_dbrec_status(dbrec_id, 'ready')
    else:
        _update_dbrec_status(dbrec_id, 'error')
    return


def get_version_before(time_stamp):
    return 'initial'


def get_genome_index_dir(gtf_ver):
    """
    Get output genome index directory for gtf_ver.
    Generated once per .gtf file.
    """

    """Return path to file or None if file not available."""
    # query DB
    mappings = enCount.db.gtfs.find({'gtf_ver': gtf_ver})

    assert(mappings.count() <= 1)
    mappings = list(mappings)

    if mappings:
        # fetch record from DB
        mapping = mappings[0]
        if mapping['status'] == 'ready':
            # genome index exists
            return mapping['out_genome_dir']
        else:
            # not ready
            return
    else:
        # GTF is in format time_stamp.gtf
        in_gtf = os.path.join(genomes_root, "gtf", gtf_ver + ".gtf")
        if not os.path.exists(in_gtf):
            print("Error, {:s} does not exist".format(in_gtf))
            return

        # Output directory
        rel_folder = os.path.join("index", gtf_ver)
        abs_folder = os.path.join(genomes_root, rel_folder)
        if not os.path.isdir(abs_folder):
            try:
                os.makedirs(abs_folder)
            except OSError:
                print('Error, could not create genome index folder: '
                      '{:s}'.format(abs_folder))
                print(' genome index {:s} will not be '
                      'generated.'.format(abs_folder))
                # error, will not be ready
                return

        time_stamp = datetime.datetime.utcnow()
        new_rec = {'gtf_ver': gtf_ver, 'status': 'to generate',
                   'time_stamp': time_stamp,
                   'out_genome_dir': abs_folder,
                   'in_gtf': in_gtf,
        }
        print('adding new record to gtfs collection: {:s}'.format(
            str(new_rec)))
        enCount.db.gtfs.insert_one(new_rec)
        # not ready yet
        return abs_folder
    return


def process(mock=False):
    """Synchronizes database and queue of current mappings jobs."""
    global submitted_gtf_generates

    # query DB to get all records that have status 'to generate'
    for e in enCount.db.gtfs.find({'status': 'to generate'}):
        genome_index_id = str(e['_id'])
        out_genome_dir = e['out_genome_dir']
        in_gtf = e['in_gtf']

        # enqueue new gtf generates to process
        # TODO: Mocking does not apply here because the method is run in a different container
        if genome_index_id not in submitted_gtf_generates:
            job = enCount.queues.gtfs.enqueue_call(
                enCount.gtfs.generate_genome_index,
                args=(in_gtf, enCount.config.genome_fasta_root, out_genome_dir, mock),
                result_ttl=-1, ttl=-1, timeout=-1,
            )
            job.meta['genome_index_id'] = genome_index_id
            job.save()
            submitted_gtf_generates[genome_index_id] = job

    # clean queue for finished mappings
    for job in enCount.queues.gtfs.jobs:
        if job.is_finished or job.is_failed:
            job.cleanup()
            enCount.queues.gtfs.remove(job)