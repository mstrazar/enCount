"""
Interface to query and download from ENCODE.

"""

import os
import gzip
import csv
import requests
import json

HEADERS = {'accept': 'application/json'}


def get_list(assay_titles=('RNA-seq', 'shRNA/RNA-seq'),
             organisms=('Homo sapiens', ),
             biosamples=('K562', 'HepG2')):

    # get json with link to batch download file
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

    print("Getting info on matching experiments...")
    url = 'https://www.encodeproject.org/search/'
    response = requests.get(url, params=params, headers=HEADERS)
    print(response.url)
    resp_json_dict = response.json()

    # get link for batch download
    batch_download_url = resp_json_dict['batch_download']

    # get link to metadata file
    response = requests.get(batch_download_url)
    metadata_url = response.text.splitlines()[0]
    host, query = metadata_url.split('/metadata/')
    query = requests.utils.quote(query)
    metadata_url = '{:s}/metadata/{:s}'.format(host, query)

    # download metadata file
    print("Downloading metadata file...")
    response = requests.get(metadata_url)
    print(response.url)
    reader = csv.reader(response.text.splitlines(), delimiter='\t')
    metadata_records = [row for row in reader]
    header = metadata_records[0]
    exps = {}
    for rec in metadata_records[1:]:
        cur = dict(zip(header, rec))
        exp_acc = cur['Experiment accession']
        exps.setdefault(exp_acc, []).append(cur)

    print('There are {0:d} experiments.'.format(len(exps)))
    return exps, metadata_records


def test():
    if os.path.isfile('tmp.json.gz'):
        resp = json.loads(gzip.open('tmp.json.gz', 'rt').read())
    else:
        resp = _get()
        gzip.open('tmp.json', 'wt').write(json.dumps(resp))
    return parse_json(resp), resp
