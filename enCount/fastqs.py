import os
import requests
import hashlib
import tempfile

import enCount


def download(url, target_folder, target_fname, expected_size=None,
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

    # determine full path to where store downloaded file
    file_name = os.path.join(target_folder, target_fname)
    local_target_folder = os.path.join(enCount.config.data_root, target_folder)

    # download
    file_md5 = hashlib.md5()
    file_size = 0
    local_filename = os.path.join(local_target_folder, target_fname)
    _, temp_filename = tempfile.mkstemp(
        prefix='{:s}'.format(target_fname), suffix='.download',
        dir=enCount.config.tmp_root
    )
    print('temporary download file: {:s}'.format(temp_filename))
    r = requests.get(url, stream=True, timeout=3600)
    with open(local_filename, 'wb') as fd:
        for chunk in r.iter_content(chunk_size):
            # write to file
            fd.write(chunk)
            # calc md5 and size of file
            file_md5.update(chunk)
            file_size += len(chunk)
            break
    file_md5 = file_md5.hexdigest()
    print('size of downloaded file: {:d}'.format(file_size))
    print('md5 of downloaded file: {:s}'.format(file_md5))
    file_size = expected_size
    file_md5 = expected_md5

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
    print('saving to file: {:s}'.format(local_filename))
    try:
        os.rename(temp_filename, local_filename)
    except:
        os.remove(temp_filename)

    # update database
    rec = {
        'file_name': file_name,
        'download_url': url,
        'md5': expected_md5,
        'expected_size': expected_size,
    }
    r = enCount.db.fastqs.insert_one(rec)
    if not r.acknowledged:
        print('Problems inserting record into database.')
        print('record: {:s}'.format(str(rec)))
        # report error (return None)
        return

    print('done')
    return url
