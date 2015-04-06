#!/usr/bin/env python

import os
import gzip
import requests

from cref.structure.secondary import SecondaryStructureBD


def download_ss(filename):
    chunk_size = 1024
    if not os.path.isfile(filename):
        response = requests.get('http://www.rcsb.org/pdb/files/ss.txt.gz')
        with open(filename, 'wb') as fd:
            for chunk in response.iter_content(chunk_size):
                fd.write(chunk)


def save_ss_to_db(filename):
    with gzip.open(filename, 'rb') as ss_file:
        i = 0
        sequence = ''
        structure = ''
        line = ss_file.readline().decode('utf-8')

        while line:
            if line.startswith('>'):
                parts = line.split(':')
                pdb = parts[0][1:]
                chain = parts[1]

                line = ss_file.readline().decode('utf-8')
                while line and not line.startswith('>'):
                    sequence += line
                    line = ss_file.readline().decode('utf-8')

                line = ss_file.readline().decode('utf-8')
                while line and not line.startswith('>'):
                    structure += line
                    line = ss_file.readline().decode('utf-8')

                structure = structure.replace('\n', '')
                structure = structure.replace(' ', '-')
                sequence = sequence.replace('\n', '')
                db.save(pdb, chain, sequence, structure)
                i += 1
                print('Saved', i, pdb, chain)
                sequence = ''
                structure = ''


if __name__ == '__main__':
    filename = "ss.txt.gz"
    db = SecondaryStructureBD()
    db.create()
    download_ss(filename)
    save_ss_to_db(filename)
    db.close()
