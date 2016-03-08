import redis
import rq

from enCount import config

# connect to Redis, needed by the queueing system rq
print('Connecting to Redis and rq queuing system host "{:s}", port: '
      '{:d}, db: {:d}'.format(config.REDIS_HOSTNAME, config.REDIS_PORT,
                              config.REDIS_DB))

_redis_conn = redis.Redis(host=config.REDIS_HOSTNAME,
                                port=config.REDIS_PORT, db=config.REDIS_DB)
# redis_conn.flushall() # call this to empty the redis database

downloads = rq.Queue('download', connection=_redis_conn, default_timeout=-1)
mappings = rq.Queue('mappings', connection=_redis_conn, default_timeout=-1)
failed = rq.get_failed_queue(connection=_redis_conn)

print(' jobs in downloads queue: {:d}'.format(downloads.count))
print(' jobs in mappings queue: {:d}'.format(mappings.count))
print(' jobs in failed queue: {:d}'.format(failed.count))
