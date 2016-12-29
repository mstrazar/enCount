import os
import importlib.machinery

from setuptools import setup, find_packages
from subprocess import check_output, CalledProcessError

here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, 'README.md')) as f:
    README = f.read()

requires = [
    'requests',
    ]


classifiers = [
    'Programming Language :: Python',
    ]


# get version from iCount/_version.py
try:
    loader = importlib.machinery.SourceFileLoader('enCount._version',
                                                  'enCount/_version.py')
    version = loader.load_module()
    version = version.__version__
except FileNotFoundError as e:
    print("Please, add a proper version file enCount/_version.py.")
    exit(1)

setup(
    name='enCount',
    version=version,
    description='Computational pipeline for analysis of ENCODE data',
    long_description=README,
    classifiers=classifiers,
    author='University of Ljubljana, Bioinformatics Laboratory',
    author_email='tomaz.curk@fri.uni-lj.si',
    url='',
    keywords='RNA-Seq iCLIP ENCODE',
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    install_requires=requires,
    tests_require=requires,
)
