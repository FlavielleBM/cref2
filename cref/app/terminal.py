#!/usr/bin/env python

import os
import argparse
import logging
import importlib
import tempfile
import subprocess

import pandas
from Bio import SeqIO

from cref.app import BaseApp
from cref.libs import rcsb

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


def parse_args():
    parser = argparse.ArgumentParser(
        description='CReF: Protein structure prediction')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '--sequence', dest='sequence',
        help='Aminoacid sequence using one letter code',
    )
    group.add_argument(
        '--fasta', dest='fasta',
        help='File containing the fasta sequence',
    )
    group.add_argument(
        '--pdb', dest='pdb',
        help='PDB Code from where the sequence will be extracted',
    )
    parser.add_argument(
        '--config', dest='config',
        help='File specifying the configurations'
    )
    parser.add_argument(
        '--output', dest='output_dir',
        default='predictions/tmp',
        help='Directory to save the results'
    )
    parser.add_argument(
        '--log', dest='log_level',
        default='INFO',
        help='Log level to be used (DEBUG, INFO, WARN, ERROR)'
    )
    parser.add_argument(
        '--pymol', dest='pymol', action='store_true',
        help='View prediction in PyMOL'
    )
    return parser.parse_args()


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
        output_dir = os.path.join(output_dir, sequence.id.split(':')[0] + '/')
        output = run_cref(seq, output_dir, params)
        sequence_file = os.path.join(output_dir, 'sequence.txt')
        with open(sequence_file, 'w') as sequence_output:
            sequence_output.write(seq)
        output_filepaths.append(output)
    return output_filepaths

def read_config(module):
    try:
        config = importlib.import_module(module)
    except Exception as e:
        logger.error(e)
        raise Exception('Invalid config file')
    return config


def run_pymol(pdb_code, predicted_filepath):
    filepath = os.path.join(
        os.path.dirname(predicted_filepath),
        'experimental_structure.pdb'
    )
    experimental_pdb = rcsb.download_pdb(pdb_code, filepath)
    subprocess.call([
        'pymol',
        predicted_filepath,
        experimental_pdb,
        '-r',
        'cref/utils/pymol.py'
    ])


def main():
    params = {}
    args = parse_args()
    configure_logger(args.log_level)
    if args.config:
        config = read_config(args.config)
        params = config.params
    # Sequence input
    if args.sequence:
        run_cref(args.sequence, args.output_dir, params)
    # Fasta file input
    elif args.fasta:
        predict_fasta(args.fasta, args.output_dir, params)
    # PDB code input
    elif args.pdb:
        handler, fasta_file = tempfile.mkstemp(suffix='.fasta', prefix='tmp')
        rcsb.download_fasta(args.pdb, fasta_file)
        params['pdb'] = args.pdb
        output_files = predict_fasta(fasta_file, args.output_dir, params)
        os.remove(fasta_file)
        if args.pymol:
            run_pymol(args.pdb, output_files[0])
    else:
        raise ValueError('You must specify a sequence, fasta file or pdb code')


if __name__ == '__main__':
    main()
