#!/usr/bin/env python

import itertools
import os
import logging
import shutil
import subprocess
import sys
import tempfile

import pandas
import requests
from Bio import SeqIO

from cref.app import BaseApp

logger = logging.getLogger('CReF')


class TerminalApp(BaseApp):
    """
    App to be run on the terminal
    """

    def reporter(self, state):
        pass


def run_cref(aa_sequence, output_dir, params):
    pandas.set_option('display.max_columns', 0)
    pandas.set_option('display.max_rows', 5)

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    app = TerminalApp(params)
    return app.run(aa_sequence, output_dir)


def configure_logger(log_level='INFO', include_pathname=False):
    logger = logging.getLogger('CReF')
    level = getattr(logging, log_level.upper(), None)
    if not isinstance(level, int):
        raise ValueError('Invalid log level: %s' % log_level)

    logger.propagate = False
    logger = logging.getLogger('CReF')
    logger.setLevel(level)

    ch = logging.StreamHandler()
    ch.setLevel(level)
    if include_pathname:
        template = ('%(asctime)s - %(name)s - %(levelname)s'
                    '(%(pathname)s, %(lineno)d)- %(message)s')
    else:
        template = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

    formatter = logging.Formatter(template, datefmt='%d/%m/%Y %I:%M:%S %p')
    ch.setFormatter(formatter)
    logger.addHandler(ch)


def read_fasta(filepath):
    records = []
    with open(filepath, 'rU') as fasta_file:
        records = list(SeqIO.parse(fasta_file, 'fasta'))
    return records


def predict_fasta(filepath, output_dir, params):
    sequences = read_fasta(filepath)
    output_filepaths = []
    for sequence in sequences:
        seq = str(sequence.seq).replace('X', '')
        output = run_cref(seq, output_dir, params)
        sequence_file = os.path.join(output_dir, 'sequence.txt')
        with open(sequence_file, 'w') as sequence_output:
            sequence_output.write(seq)
        output_filepaths.append(output)
    return output_filepaths


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
    url = ('http://www.rcsb.org/pdb'
           '/files/fasta.txt?structureIdList=' + pdb_code.upper())
    return _download_file(url, filepath)


def download_pdb(pdb_code, filepath):
    url = ('http://www.rcsb.org/pdb/download/downloadFile.do?'
           'fileFormat=pdb&compression=NO&structureId=' + pdb_code.upper())
    return _download_file(url, filepath)


def read_config():
    pdb_id = sys.argv[1].upper()
    excluded_pdbs = [x.strip() for x in sys.argv[2].split()]
    excluded_pdbs.append(pdb_id)

    fragment_size = [5, 7, 9]
    number_of_clusters = [4, 6, 8]
    matrix = ["PAM30", "BLOSUM62"]
    max_templates = [50, 100]
    number_of_alignments = [500]
    params_list = []
    print(len(list(itertools.product(
            fragment_size,
            number_of_clusters,
            matrix,
            max_templates,
            number_of_alignments))))
    for f, c, m, t, a in itertools.product(
            fragment_size,
            number_of_clusters,
            matrix,
            max_templates,
            number_of_alignments):
        print(f, c, t, a, m)
        params = {
            "id": '_'.join([str(x) for x in (f, c, t, a, m)]),
            "exclude": {"pdbs": excluded_pdbs},
            "fragment_size": f,
            "number_of_clusters": c,
            "max_templates": t,
            "blast": {
                "number_of_alignments": a,
                "scoring": {
                    "matrix": m,
                    "gap_costs": "ungapped",
                }
            }
        }
        params['pdb'] = pdb_id
        params['output_dir'] = os.path.join(
                'predictions/benchmark',
                params['pdb'],
                '_'.join([str(x) for x in (f, c, t, a, m)]),
        )
        params_list.append(params)
    return params_list


def run_pymol(pdb_code, predicted_filepath):
    filepath = os.path.join(
        os.path.dirname(predicted_filepath),
        'experimental_structure.pdb'
    )
    experimental_pdb = download_pdb(pdb_code, filepath)
    output = subprocess.check_output([
        'pymol',
        predicted_filepath,
        experimental_pdb,
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
    configure_logger('WARN')
    test_cases = read_config()
    results = {}

    for params in test_cases:
        try:
            print('Predicting', params['pdb'], params['id'])
            handler, fasta_file = tempfile.mkstemp(
                    suffix='.fasta', prefix='tmp')
            download_fasta(params['pdb'], fasta_file)
            output_files = predict_fasta(
                    fasta_file, params['output_dir'], params)
            rmsd, imagepath = run_pymol(params['pdb'], output_files[0])
            output_file = os.path.join(params['output_dir'], 'rmsd.txt')
            with open(output_file, 'w') as rmsd_file:
                rmsd_file.write(rmsd)
            results[params['id']] = rmsd

            shutil.copyfile(
                    imagepath,
                    os.path.join(params['output_dir'], 'alignment-pymol.png'),
            )
            print('Prediction written to', output_files)
            os.remove(fasta_file)
        except Exception as error:
            print(error)
            print(params)
    print(results)
    rmsds_to_csv(results, params['pdb'])


if __name__ == '__main__':
    main()
