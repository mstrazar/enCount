import os


# folder configuration
_config_root = os.path.split(os.path.abspath(__file__))[0]

volume_root = "/endata/"
data_root = os.path.join(volume_root, "data")
genomes_root = os.path.join(volume_root, "genomes")
results_root = os.path.join(volume_root, "results")

mappings_root = os.path.join(results_root, "mappings")
counts_root = os.path.join(results_root, "counts")

tmp_root = '/tmp/enCount'
data_debug_root = os.path.join(tmp_root, 'debug')

# Redis store
REDIS_HOSTNAME = 'redis'
REDIS_PORT = 6379
REDIS_DB = 0

# MongoDB
MONGO_HOSTNAME = 'mongodb://mongo'
MONGO_PORT = 27017

# QoRTS and JunctionSeq scripts
QORTS_JAR = os.path.join(_config_root, "externals", "libs", "QoRTs.jar")
QORTS_R = os.path.join(_config_root, "externals", "libs", "QoRTs.R")
JUNCTIONSEQ_R = os.path.join(_config_root, "externals", "libs", "JunctionSeq.R")

# STAR aligner
RSCRIPT = "/usr/bin/Rscript"
STAR_EXEC = "/home/enuser/bin/STAR"
NUM_THREADS = 4
RAM_LIMIT = None

# DEXSEQ requirements
PYTHON2_EXEC = "/usr/bin/python"
DEXSEQ_PREP_ANNOTATION = "/home/enuser/.R/DEXSeq/python_scripts/dexseq_prepare_annotation.py"

# Test data locations
ENDATA_TEST_REPO = "https://github.com/mstrazar/enCount-data"


# read (override) from local config
try:
    from .myconfig import *
except ImportError:
    print('enCount: You can define a local configuration in myconfig.py, '
          'which needs to be stored in folder: {:s}.'.format(_config_root))