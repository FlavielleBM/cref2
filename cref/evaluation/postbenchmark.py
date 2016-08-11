#!/usr/bin/env python

import itertools
import os
import shutil
import subprocess
import sys

import pandas
import requests


def _download_file(url, filepath):
    r = requests.get(url, stream=True)
    with open(filepath, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024):
            if chunk:  # filter out keep-alive new chunks
                f.write(chunk)
                f.flush()
    return filepath


def download_pdb(pdb_code, filepath):
    url = ('http://www.rcsb.org/pdb/download/downloadFile.do?'
           'fileFormat=pdb&compression=NO&structureId=' + pdb_code.upper())
    return _download_file(url, filepath)


def read_config():
    pdb_id = sys.argv[1].upper()
    fragment_size = [5, 7, 9]
    number_of_clusters = [4, 6, 8, 10]
    matrix = ["PAM30", "BLOSUM62"]
    max_templates = [50, 100, 200]
    number_of_alignments = [500, 1000]
    params_list = []
    for f, c, m, t, a in itertools.product(
            fragment_size,
            number_of_clusters,
            matrix,
            max_templates,
            number_of_alignments):
        params = {}
        params['id'] = (f, c, t, a, m)
        params['pdb'] = pdb_id
        params['output_dir'] = os.path.join(
                'predictions/benchmark',
                params['pdb'],
                '_'.join([str(x) for x in (f, c, t, a, m)]),
        )
        params_list.append(params)
    return params_list


def run_pymol(pdb_code, output_dir):
    predicted_filepath = os.path.join(
        output_dir, 'predicted_structure.pdb')
    experimental_filepath = os.path.join(
        output_dir, 'experimental_structure.pdb')
    if not os.path.isfile(experimental_filepath):
        download_pdb(pdb_code, experimental_filepath)

    output = subprocess.check_output([
        'pymol',
        predicted_filepath,
        experimental_filepath,
        '-r',
        'cref/utils/pymolbench.py'
    ])
    output = output.decode('utf-8').split('\n')
    rmsd = output[-4]
    imagepath = output[-3]
    return rmsd, imagepath


def rmsds_to_csv(rmsds, filename):
    results = []

    for key, value in rmsds.items():
            results.append(key + (value,))

    df = pandas.DataFrame(
        results,
        columns=['fragment_size', 'group_count', 'max_templates',
                 'max_blast', 'matrix', 'rmsd']
    )
    df.to_csv(filename + '.rmsd.csv')


def main():
    test_cases = read_config()
    results = {}

    for params in test_cases:
        try:
            rmsd, imagepath = run_pymol(params['pdb'], params['output_dir'])
            results[params['id']] = rmsd

            output_file = os.path.join(params['output_dir'], 'rmsd.txt')
            with open(output_file, 'w') as rmsd_file:
                rmsd_file.write(rmsd)

            shutil.copyfile(
                    imagepath,
                    os.path.join(params['output_dir'], 'alignment-pymol.png'),
            )
        except Exception as error:
            print(error)
            print(params)

    print(results)
    rmsds_to_csv(results, params['pdb'])

if __name__ == '__main__':
    main()
