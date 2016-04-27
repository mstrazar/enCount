"""
Interface to query and download from ENCODE.

"""

import csv
import requests
import urllib

HEADERS = {'accept': 'application/json'}


def parse_metadata_records(metadata_records):
    """Group metadata records by Experiment accession"""
    header = metadata_records[0]
    exps = {}
    for rec in metadata_records[1:]:
        cur = dict(zip(header, rec))
        exp_acc = cur['Experiment accession']
        exps.setdefault(exp_acc, []).append(cur)

    print('There are {0:d} experiments.'.format(len(exps)))
    return exps


def get_online_list(assay_titles=('RNA-seq', 'shRNA+RNA-seq'),
                    organisms=('Homo sapiens', ),
                    biosamples=('K562', 'HepG2')):

    params = {
        'type': 'Experiment',
        'assay_title': assay_titles,
        'biosample_term_name': biosamples,
        'replicates.library.biosample.donor.organism.scientific_name':
            organisms,
        'status': 'released',
        'files.file_type': 'fastq',
        'limit': 'all',
        'frame': 'objects',
    }

    print('Getting info on matching experiments...')
    params = urllib.parse.urlencode(params, doseq=True)
    params = urllib.parse.quote(params)
    url = 'https://www.encodeproject.org/metadata/{:s}/metadata.csv'.format(
        params)
    response = requests.get(url)
    print(response.url)
    reader = csv.reader(response.text.splitlines(), delimiter='\t')
    metadata_records = [row for row in reader]
    print('Retrieved {:d} records.'.format(len(metadata_records)))
    return parse_metadata_records(metadata_records)
