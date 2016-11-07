import os
import enCount
import datetime
import time
import hashlib
from bson.objectid import ObjectId
from enCount.config import genomes_root

submitted_gtf_generates = dict(
    (j.meta['genome_index_id'], j) for j in enCount.queues.gtfs.jobs
)

def _update_dbrec_status(dbrec_id, new_status,):
    """
    Update a database record on a .gtf file.
    """
    r = enCount.db.mappings.update_one(
        {'_id': ObjectId(dbrec_id)}, {"$set": {'status': new_status,
            'out_genome_dir': out_genome_dir}}
    )
    if not r.acknowledged:
        print(' problems updating collection mappings record id: {'
              ':s}'.format(dbrec_id))


def generate_genome_index(in_gtf, in_genome_fasta_dir, dbrec_id):
    """
    Generate a genome index using STAR for a given .gtf file.
    """
    gtf_ver = get_version_before(datetime.datetime.utcnow())
    out_genome_dir = enCount.gtfs.get_genome_index_dir(gtf_ver=gtf_ver)

    r = enCount.externals.rnastar.run_star_generate_genome(in_gtf=in_gtf,
                                 in_genome_fasta_dir=in_genome_fasta_dir,
                                 out_genome_dir=out_genome_dir,
                                 num_threads=enCount.config.NUM_THREADS)
    if r == 0:
        _update_dbrec_status(dbrec_id, 'ready', gtf_ver)
    else:
        _update_dbrec_status(dbrec_id, 'error', gtf_ver)
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
            return mapping['']
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
            except:
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
        # not ready
        return abs_folder
    return



def process():
    """Synchronizes database and queue of current mappings jobs."""
    global submitted_gtf_generates
    # query DB to get all records that have status 'to generate'
    for e in enCount.db.gtfs.find({'status': 'to generate'}):
        genome_index_id = str(e['_id'])
        gtf_ver = e['gtf_ver']
        out_genome_dir = e['out_genome_dir']
        in_gtf = e['in_gtf']


        # enqueue new gtf generates to process
        if genome_index_id not in submitted_gtf_generates:

            job = enCount.queues.gtfs.enqueue_call(
                enCount.gtfs.generate_genome_index,
                args=(in_gtf, enCount.config.genome_fasta_root, out_genome_dir),
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

if __name__ == "__main__":
    # Get initial, minimal version
    gtf_ver = get_version_before(datetime.datetime.min)
    assert os.path.exists(os.path.join(genomes_root, "gtf", "%s.gtf" % gtf_ver))

    # Create new genome index directory
    get_genome_index_dir(gtf_ver)
    assert os.path.isdir(os.path.join(genomes_root, "index", "initial"))

    print("Initial queue status")
    enCount.queues.gtfs.empty()
    enCount.queues.print_stats()

    # Run process loop
    process()

    is_finished = False
    while not is_finished:
        time.sleep(5)
        print("Waiting for genome index generation to finish")
        enCount.queues.print_stats()
        is_finished = enCount.queues.gtfs.is_empty()
        print()

    # Assert results exist
    outdir = get_genome_index_dir(gtf_ver)
    assert outdir is not None
    assert os.path.exists(os.path.join(outdir, "genomeParameters.txt"))