"""
Interface to query and download from ENCODE.

"""

import csv
import requests
import urllib

HEADERS = {'accept': 'application/json'}


def parse_metadata_records(metadata_records, header):
    """ Group metadata records by Experiment accession

        metadata_records is a dict indexed by fastq_ID.

        If experiment is a knockdown, add related control experiments.
        There may be multiple knockdown experiments with same controls, so controls are
        stored with a different experiment ID.

        'Controlled by' filed holds pointer to files, e.g.
        '/files/ENCFF078MXU/, /files/ENCFF791HTS/'
        that point to fastq IDs directly (and not experiment IDs).

        Important: Make sure experiment IDs are not used downstream of here because they will remain stored in
        these rows.

    """
    exps = {}

    for _, rec in metadata_records.items():
        # Group fastq. file entries by their experiment ID
        cur = dict(zip(header, rec))
        exp_acc = cur['Experiment accession']
        exps.setdefault(exp_acc, []).append(cur)

        # Controls
        if cur['Controlled by'] != "":
            control_ids = map(lambda e: e.split("/")[2], cur['Controlled by'].split(","))
            for cid in control_ids:
                # Find rows for these files
                if cid in metadata_records:
                    crow = dict(zip(header, metadata_records[cid]))
                    if crow not in exps[exp_acc]:
                        exps[exp_acc].append(crow)
                else:
                    raise ValueError("Controls for experiment %s recorded, but not fetched!" % exp_acc)

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
    header = next(reader)
    metadata_records = dict(((row[0], row) for row in reader))
    print('Retrieved {:d} records.'.format(len(metadata_records)))
    return parse_metadata_records(metadata_records, header)
