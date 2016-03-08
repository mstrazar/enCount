import pymongo

from enCount import config

# mongoDB for storing results of processing
print('Connecting to MongoDB host "{:s}", '
      'port: {:d}'.format(config.MONGO_HOSTNAME, config.MONGO_PORT))

_client = pymongo.MongoClient(config.MONGO_HOSTNAME, config.MONGO_PORT)
db_encode = _client['encode']

fastqs = db_encode['fastqs']  # info on downloaded files
fastqs.create_index([('file_name', pymongo.ASCENDING)], unique=True)

experiments = db_encode['experiments']  # info on new metadata versions
experiments.create_index([('time_stamp', pymongo.DESCENDING)], unique=True)

mappings = db_encode['mappings']  # info on mappings

# reset collections
# col_fastqs.remove({})
# col_experiments.remove({})
# col_mappings.remove({})
print(' records on fastq files in DB: {:d}'.format(fastqs.find().count()))
print(' records on experiments in DB: {:d}'.format(experiments.find().count()))
print(' records on mappings in DB: {:d}'.format(mappings.find().count()))