#!/usr/bin/env python
import sys
import logging
from math import floor
from cref import sequence
from cref.sequence.alignment import Blast
from cref.structure import PDB, torsions


def run_cref(aa_sequence):
    pdb_downloader = PDB.PDBDownloader('data/pdb')
    blast = Blast(db='tests/blastdb/pdbseqres')
    size = 5
    central_aa = floor(size / 2)
    for fragment in sequence.fragment(aa_sequence, size):
        blast_results = list(blast.align(fragment))[:30]
        aa = []
        phi = []
        psi = []
        for result in blast_results:
            try:
                pdb_downloader.retrieve(result.pdb_code)
                print(result.pdb_code)
                print(result.chain)
                torsion_angles = torsions.backbone_torsion_angles(
                    result.pdb_code,
                    result.chain,
                    'data/pdb/pdb{}.ent'.format(result.pdb_code)
                )
                for hit in result.hits:
                    target_aa = central_aa - hit['query_start'] + 1
                    subject_start = torsion_angles['residues'].find(
                        hit['subject'])
                    if target_aa >= 0 and subject_start >= 0 and \
                            hit['subject'][central_aa] == hit['query'][central_aa]:
                        target_aa = subject_start + target_aa
                        phi.append(torsion_angles['phi'][target_aa])
                        psi.append(torsion_angles['psi'][target_aa])
                        aa.append(torsion_angles['residues'][target_aa])
            except Exception as error:
                logging.error("Could not download " + result.pdb_code)
                logging.error(error)
        print(aa, phi, psi)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        run_cref(sys.argv[1])
    else:
        print('Syntax: cref <aminoacid sequence>')
