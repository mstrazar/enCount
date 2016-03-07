"""
enCount tasks and analyses.

enCount is a Python library for processing RNA-Seq data from ENCODE.

"""


import os

# from ._version import __version__
from . import encode
from . import externals
from . import workers

# CONFIG
# storage folder location
_config_root = os.path.split(os.path.abspath(__file__))[0]

data_root = os.path.join(_config_root, 'data')
genomes_root = os.path.join(_config_root, 'genomes')
results_root = os.path.join(_config_root, 'results')
tmp_root = os.path.join(_config_root, 'tmp')

# read from local config
try:
    from .myconfig import *
except ImportError:
    print('enCount: Warning, please define local config in myconfig.py and '
          'store it to folder: {:s}.'.format(_config_root))

# create folders if not present yet
if not os.path.exists(data_root):
    print('Data root folder does not exist. Will create it at: '
          '{:s}'.format(data_root))
    os.makedirs(data_root)

if not os.path.exists(genomes_root):
    print('Genomes root folder does not exist. Will create it at: '
          '{:s}'.format(genomes_root))
    os.makedirs(genomes_root)

if not os.path.exists(results_root):
    print('Results root folder does not exist. Will create it at: '
          '{:s}'.format(results_root))
    os.makedirs(results_root)

if not os.path.exists(tmp_root):
    print('Temporary files and folders root folder does not exist. '
          'Will create it at: {:s}'.format(tmp_root))
    os.makedirs(tmp_root)
