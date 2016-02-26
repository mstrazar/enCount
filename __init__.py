import os

_config_root = os.path.split(os.path.abspath(__file__))[0]

data_root = os.path.join(_config_root, 'data')
genomes_root = os.path.join(_config_root, 'genomes')
results_root = os.path.join(_config_root, 'results')

# read from local config
try:
    from .myconfig import *
except ImportError:
    print(u'enCount: Warning, please define local config in myconfig.py and '
          u'store it to folder: {0:s}.'.format(_config_root))
