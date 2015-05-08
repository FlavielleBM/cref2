#!/usr/bin/env python

import os
import argparse
import logging
import importlib
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


def run_cref(aa_sequence, output_dir, params):
    pandas.set_option('display.max_columns', 0)
    pandas.set_option('display.max_rows', 5)

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    app = TerminalApp(params)
    return app.run(aa_sequence, output_dir)


def configure_logger(log_level='INFO'):
    logger = logging.getLogger('CReF')
    level = getattr(logging, log_level.upper(), None)
    if not isinstance(level, int):
        raise ValueError('Invalid log level: %s' % log_level)

    logger.propagate = False
    logger = logging.getLogger('CReF')
    logger.setLevel(level)

    ch = logging.StreamHandler()
    ch.setLevel(level)
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%d/%m/%Y %I:%M:%S %p'
    )

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
        '--output_dir', dest='output_dir',
        default='predictions/tmp',
        help='Directory to save the results'
    )
    parser.add_argument(
        '--log', dest='log_level',
        default='INFO',
        help='Log level to be used (DEBUG, INFO, WARN, ERROR)'
    )
    return parser.parse_args()


def read_fasta(filepath):
    records = []
    with open(filepath, 'rU') as fasta_file:
        records = list(SeqIO.parse(fasta_file, 'fasta'))
    return records


def predict_fasta(filepath, output_dir, params):
    sequences = read_fasta(filepath)
    for seq in sequences:
        # Remove unknown aminoacids
        sequence = str(seq.seq).replace('X', '')
        run_cref(
            sequence,
            os.path.join(output_dir, seq.id.split('|')[0]),
            params
        )


def download_fasta(pdb_code, filepath):
    """"""
    url = ('http://www.rcsb.org/pdb'
           '/files/fasta.txt?structureIdList=' + pdb_code.upper())
    r = requests.get(url, stream=True)
    with open(filepath, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024):
            if chunk:  # filter out keep-alive new chunks
                f.write(chunk)
                f.flush()
    return filepath


def read_config(module):
    try:
        config = importlib.import_module(module)
    except Exception as e:
        logger.error(e)
        raise Exception('Invalid config file')
    return config


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
        download_fasta(args.pdb, fasta_file)
        predict_fasta(fasta_file, args.output_dir, params)
        os.remove(fasta_file)
    else:
        raise ValueError('You must specify a sequence, fasta file or pdb code')


if __name__ == '__main__':
    main()
