import requests

def _download_file(url, filepath):
    r = requests.get(url, stream=True)
    with open(filepath, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024):
            if chunk:  # filter out keep-alive new chunks
                f.write(chunk)
                f.flush()
    return filepath


def download_fasta(pdb_code, filepath):
    """"""
    url = ('https://www.rcsb.org/pdb/download/downloadFastaFiles.do?'
           'structureIdList={}&compressionType=uncompressed'.format(pdb_code.upper()))
    return _download_file(url, filepath)


def download_pdb(pdb_code, filepath):
    url = 'https://files.rcsb.org/download/{}.pdb'.format(pdb_code.upper())
    return _download_file(url, filepath)
