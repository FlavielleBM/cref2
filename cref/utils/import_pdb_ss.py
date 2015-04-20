#!/usr/bin/env python

import os
import requests

from cref.structure.secondary import PDBSecondaryStructureDB


def download_ss(filename):
    chunk_size = 1024
    if not os.path.isfile(filename):
        response = requests.get('http://www.rcsb.org/pdb/files/ss.txt.gz')
        with open(filename, 'wb') as fd:
            for chunk in response.iter_content(chunk_size):
                fd.write(chunk)


def save_ss_to_db(filename):
    with open(filename, 'r') as ss_file:
        i = 0
        sequence = ''
        structure = ''
        line = ss_file.readline()

        while line:
            if line.startswith('>'):
                parts = line.split(':')
                pdb = parts[0][1:]
                chain = parts[1]

                line = ss_file.readline()
                while line and not line.startswith('>'):
                    sequence += line
                    line = ss_file.readline()

                line = ss_file.readline()
                while line and not line.startswith('>'):
                    structure += line
                    line = ss_file.readline()

                structure = structure.replace('\n', '')
                structure = structure.replace(' ', '-')
                sequence = sequence.replace('\n', '')
                db.save(pdb, chain, sequence, structure)
                i += 1
                print('Saved', i, pdb, chain)
                sequence = ''
                structure = ''


if __name__ == '__main__':
    import sys
    db = PDBSecondaryStructureDB(sys.argv[2])
    db.create()
    save_ss_to_db(sys.argv[1])
    db.close()
