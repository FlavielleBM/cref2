#!/usr/bin/env python
import sys
import logging
from cref import sequence
from cref.sequence.alignment import Blast
from cref.structure import PDB, torsions


def run_cref(aa_sequence):
    pdb_downloader = PDB.PDBDownloader('data/pdb')
    blast = Blast(db='tests/blastdb/pdbseqres')
    for fragment in sequence.fragment(aa_sequence):
        print(fragment)
        blast_results = blast.align(fragment)
        for result in blast_results:
            try:
                pdb_downloader.retrieve(result.pdb_code)
                torsion_angles = torsions.backbone_torsion_angles(
                    'data/pdb/pdb{}.ent'.format(result.pdb_code))
                for hit in result.hits:
                    print(hit)
                    print(
                        torsion_angles[hit['subject_start'] - 1:
                                       hit['subject_end']])
            except Exception as error:
                logging.error("Could not download " + result.pdb_code)
                logging.error(error)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        run_cref(sys.argv[1])
    else:
        print('Syntax: cref <aminoacid sequence>')
