import os
import enCount
import hashlib
import datetime
from bson.objectid import ObjectId
import sys


# Constantly named files
STAR_BAM_NAME = "Aligned.sortedByCoord.out.bam"
STAR_BAM_LOG = "Log.final.out"
QORTS_COUNT_NAME = "QC.spliceJunctionCounts.knownSplices.txt.gz"

# populate list with currently queued jobs
submitted_mappings = dict(
    (j.meta['mapping_id'], j) for j in enCount.queues.mappings.jobs
)


def _update_dbrec_status(dbrec_id, new_status, new_bam=None, read_count=None):
    """
    Update a status in the mappings database.
    :param dbrec_id: Record ID.
    :param new_status: New status: error, to count, ready.
    :param new_bam: Optionally, new BAM file .
    :param read_count: Optionally fastq read count.
    """
    row = {'status': new_status}
    if new_status == "to count": row['out_bam_file'] = new_bam
    if new_status == "to count": row["read_count"] = read_count

    r = enCount.db.mappings.update_one({'_id': ObjectId(dbrec_id)}, {"$set": row})
    if not r.acknowledged:
        print(' problems updating collection mappings record id: {0:s}'.format(dbrec_id))



def map_fastq(fastq_pair, in_genome_dir, out_dir, dbrec_id):

    enCount.externals.rnastar.run_star(
        in_genome_dir=in_genome_dir,
        in_fastq_pair=fastq_pair,
        out_dir=out_dir,)

    new_bam = os.path.join(out_dir, STAR_BAM_NAME)
    log_final = os.path.join(out_dir, STAR_BAM_LOG)
    read_count = enCount.externals.rnastar.get_read_count(out_dir)

    if os.path.exists(log_final):
        _update_dbrec_status(dbrec_id, 'to count', new_bam, read_count)
    else:
        _update_dbrec_status(dbrec_id, 'error', new_bam, read_count)
    return


def count_bam(in_bam_dir, in_gtf, out_count_dir, dbrec_id):
    """
    Map mapped reads in BAM to counting bins in a GFF.

    :param in_bam_dir: Input BAM file.
    :param in_gff: GFF bin file.
    :param out_count_dir: Output directory.
    :return:
    """
    r = enCount.externals.junctionseq.run_QoRTs_count(in_bam=in_bam_dir,
                                                      in_gtf=in_gtf,
                                                      out_dir=out_count_dir)
    if r == 0:
        if dbrec_id: _update_dbrec_status(dbrec_id, 'ready')
    else:
        if dbrec_id: _update_dbrec_status(dbrec_id, 'error')
    return


def get_mapping_id(fastq_pair, gtf_ver):
    # create a mapping ID tied to fastq_pair and .gtf version
    base_str = "%s.%s" % (fastq_pair, gtf_ver)
    return hashlib.md5(str(base_str).encode("utf-8")).hexdigest()


def get_mapping_data(bam_path, field):
    """
    Get arbitrary field value from mapping database for ready BAM/count files.

    :param bam_path: Path to bam file.
    :param field: Database field.
    :return:
    """
    mappings = enCount.db.mappings.find({'out_dir': bam_path, 'status': "ready"})
    if mappings.count() == 1:
        return next(mappings)[field]
    else:
        return None


def get_bam_file_paths(fastq_pair, gtf_ver):

    """Return path to file or None if file not available."""
    # query DB
    mappings = enCount.db.mappings.find({'fastq_pair': fastq_pair, 'gtf_ver': gtf_ver})
    assert(mappings.count() <= 1)
    mappings = list(mappings)

    if mappings:
        # fetch record from DB
        mapping = mappings[0]
        # return path only if completely ready (mapped and counted), so queues are synced
        if mapping['status'] == 'ready':
            return mapping['out_dir']
        else:
            # not ready
            return None
    else:
        # Create a new mapping directory
        map_id   = get_mapping_id(fastq_pair, gtf_ver)
        out_dir = os.path.join(enCount.config.mappings_root, gtf_ver, map_id)
        out_count_dir = os.path.join(enCount.config.counts_root, gtf_ver, map_id)

        if not os.path.exists(out_dir): os.makedirs(out_dir)
        if not os.path.exists(out_count_dir): os.makedirs(out_count_dir)

        # Check for genome index directory
        gtf_index = enCount.gtfs.get_genome_index_dir(gtf_ver)
        if gtf_index is None:
            # not ready, not inserted
            return None

        time_stamp = datetime.datetime.utcnow()
        new_rec = {'fastq_pair': fastq_pair, 'gtf_ver': gtf_ver, 'gtf_index': gtf_index,
                   'status': 'to map', 'time_stamp': time_stamp,
                   'out_dir': out_dir, 'read_count': -1, 'out_count_dir': out_count_dir,
        }
        print('adding new record to mappings collection: {:s}'.format(
            str(new_rec)))
        enCount.db.mappings.insert_one(new_rec)
        # not ready
        return None


def get_count_file_paths(bam_dir, gtf_ver):
    """
     Get count directory. Must be called after get_bam_file_paths not to duplicate inputs

    :param bam_dir: BAM directory on output.
    :param gtf_ver: GTF version.
    :return:
    """
    mappings = enCount.db.mappings.find({'out_dir': bam_dir, 'gtf_ver': gtf_ver})
    if mappings.count():
        return next(mappings)["out_count_dir"]
    else:
        return None


def process():
    """Synchronizes database and queue of current mappings jobs."""
    global submitted_mappings

    # query DB to get all records that have status 'to map' or 'to_count'
    for e in enCount.db.mappings.find({"$or": [{'status': 'to map'}, {'status': 'to count'}]}):
        dbrec_id = str(e['_id'])
        out_bam_dir = e['out_dir']
        out_count_dir = e['out_count_dir']
        status = e["status"]
        fastq_pair = e['fastq_pair']
        gtf_index = e['gtf_index']
        gtf_ver = e['gtf_ver']
        job_id = "%s_%s" % (dbrec_id, status)

        in_gtf = enCount.gtfs.version_to_path(gtf_ver)
        in_bam = os.path.join(out_bam_dir, STAR_BAM_NAME)

        # queue new mappings to process
        if job_id not in submitted_mappings:
            if status == "to map":
                job = enCount.queues.mappings.enqueue_call(
                    enCount.mappings.map_fastq,
                    args=(fastq_pair, gtf_index, out_bam_dir, dbrec_id),
                    result_ttl=-1, ttl=-1, timeout=-1,)
            else:
                job = enCount.queues.mappings.enqueue_call(
                    enCount.mappings.count_bam,
                    args=(in_bam, in_gtf, out_count_dir, dbrec_id),
                    result_ttl=-1, ttl=-1, timeout=-1,)

            job.meta['mapping_id'] = job_id
            job.save()
            submitted_mappings[job_id] = job


    # clean queue for finished mappings
    for job in enCount.queues.mappings.jobs:
        if job.is_finished or job.is_failed:
            job.cleanup()
            enCount.queues.mappings.remove(job)
