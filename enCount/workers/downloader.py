import os
import requests
import hashlib
import tempfile

from pymongo import MongoClient


client = MongoClient()
db_encode = client['encode']
col_fastqs = db_encode['fastqs']  # info on downloaded files


def fastq_download(url, target_folder, target_fname, expected_size=None,
                   expected_md5=None, chunk_size=20000000):
    print('Downloading from: {:s}'.format(url))
    if expected_size is None:
        print('no expected size specified')
    else:
        print('expected size: {:d}'.format(expected_size))
    if expected_md5 is None:
        print('no expected md5 specified')
    else:
        print('expected md5: {:s}'.format(expected_md5))

    # download
    file_md5 = hashlib.md5()
    file_size = 0
    filename = os.path.join(target_folder, target_fname)
    _, temp_filename = tempfile.mkstemp(
        prefix='.{:s}'.format(target_fname), suffix='.download',
        dir=target_folder
    )
    print('temporary file: {:s}'.format(temp_filename))
    r = requests.get(url, stream=True, timeout=3600)
    with open(temp_filename, 'wb') as fd:
        for chunk in r.iter_content(chunk_size):
            # write to file
            fd.write(chunk)
            # calc md5 and size of file
            file_md5.update(chunk)
            file_size += len(chunk)
    file_md5 = file_md5.hexdigest()
    print('size of downloaded file: {:d}'.format(file_size))
    print('md5 of downloaded file: {:s}'.format(file_md5))

    # check for errors in size or md5
    if expected_size is not None and expected_size != file_size:
        print('ERROR, size of downloaded file not as expected.')
        os.remove(temp_filename)
        return

    if expected_md5 is not None and expected_md5 != file_md5:
        print('ERROR, md5 of downloaded file not as expected.')
        os.remove(temp_filename)
        return

    # rename file to proper name
    print('saving to file: {:s}'.format(filename))
    try:
        os.rename(temp_filename, filename)
    except:
        os.remove(temp_filename)

    # update database
    rec = {
        'file_folder': target_folder,
        'file_name': target_fname,
        'download_url': url,
        'md5': expected_md5,
        'expected_size': expected_size,
    }
    r = col_fastqs.insert_one(rec)
    if not r.acknowledged:
        print('Problems inserting record into database.')
        print('record: {:s}'.format(str(rec)))
        # report error (return None)
        return

    print('done')
    return url
