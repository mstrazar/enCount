import redis
import rq

from enCount import config

# connect to Redis, needed by the queueing system rq
print('Connecting to Redis and rq queuing system host "{:s}", port: '
      '{:d}, db: {:d}'.format(config.REDIS_HOSTNAME, config.REDIS_PORT,
                              config.REDIS_DB))

_redis_conn = redis.Redis(host=config.REDIS_HOSTNAME,
                                port=config.REDIS_PORT, db=config.REDIS_DB)
# _redis_conn.flushall() # call this to empty the redis database

downloads = rq.Queue('downloads', connection=_redis_conn, default_timeout=-1)
experiments = rq.Queue('experiments', connection=_redis_conn, default_timeout=-1)
mappings = rq.Queue('mappings', connection=_redis_conn, default_timeout=-1)
failed = rq.get_failed_queue(connection=_redis_conn)


def queue_stats(q):
    cn_size = 0
    cn_queued = 0
    cn_started = 0
    cn_finished = 0
    cn_failed = 0
    for j in q.jobs:
        cn_size += 1
        cn_queued += j.is_queued
        cn_started += j.is_started
        cn_finished += j.is_finished
        cn_failed += j.is_failed
    return 'size {:d}, queued {:d}, started {:d}, finished {:d}, ' \
           'failed {:d}'.format(cn_size, cn_queued, cn_started, cn_finished,
                                cn_failed)

def print_stats():
    print(' downloads queue   : {:s}'.format(queue_stats(downloads)))
    print(' experiments queue : {:s}'.format(queue_stats(experiments)))
    print(' mappings queue    : {:s}'.format(queue_stats(mappings)))
    print(' failed queue      : {:s}'.format(queue_stats(failed)))

print_stats()