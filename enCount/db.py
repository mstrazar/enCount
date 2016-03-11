import pymongo

from enCount import config

# mongoDB for storing results of processing
print('Connecting to MongoDB host "{:s}", '
      'port: {:d}'.format(config.MONGO_HOSTNAME, config.MONGO_PORT))

_client = pymongo.MongoClient(config.MONGO_HOSTNAME, config.MONGO_PORT)
# _client.drop_database('encode')
db_encount = _client['encode']

fastqs = db_encount['fastqs']  # info on downloaded files
fastqs.create_index([
    ('e_acc', pymongo.ASCENDING),
    ('f_acc', pymongo.ASCENDING),
    ('f_md5', pymongo.ASCENDING),
    ('f_size', pymongo.ASCENDING),
], unique=True)

gtfs = db_encount['gtf']

experiments = db_encount['experiments']  # info on new metadata versions
experiments.create_index([('time_stamp', pymongo.DESCENDING)], unique=True)

mappings = db_encount['mappings']  # info on mappings

# reset collections
# fastqs.remove({})
# gtfs.remove({})
# experiments.remove({})
# mappings.remove({})
print(' records on fastq files in DB: ready {:d}/{:d}'.format(
    fastqs.find({'status': 'ready'}).count(), fastqs.find().count())
)
print(' records on gtf files in DB: ready {:d}/{:d}'.format(
    gtfs.find({'status': 'ready'}).count(), gtfs.find().count())
)
print(' records on experiments in DB: ready {:d}/{:d}'.format(
    experiments.find({'status': 'ready'}).count(), experiments.find().count())
)
print(' records on mappings in DB: ready {:d}/{:d}'.format(
    mappings.find({'status': 'ready'}).count(), mappings.find().count())
)
