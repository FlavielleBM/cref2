#!/usr/bin/env python

import os
import sys
import logging
import warnings
from math import floor

from cref import sequence
from cref.sequence.alignment import Blast
from cref.structure import PDB, torsions, plot
from cref.structure import predict_secondary_structure, write_pdb
from cref.structure.clustering import cluster_torsion_angles


class TerminalApp:

    def __init__(self):
        self.pdb_downloader = PDB.PDBDownloader('data/pdb')
        self.blast = Blast(db='tests/blastdb/pdbseqres')
        self.fragment_size = 5
        self.central = floor(self.fragment_size / 2)

    def get_central_angles(self, angles, hit):
        central = self.central
        pos = central - hit['query_start'] + 1
        subject_start = angles['residues'].find(hit['subject'])

        if pos >= 0 and subject_start >= 0 and \
                hit['subject'][central] == hit['query'][central]:
            pos = subject_start + pos
            return (
                angles['residues'][pos], angles['psi'][pos], angles['phi'][pos])
        return ('', 0, 0)

    def run(self, aa_sequence):
        fragment_angles = []
        for fragment in sequence.fragment(aa_sequence, self.fragment_size):
            torsion_angles = dict(residues='', phi=[], psi=[])
            blast_results = self.blast.align(fragment)

            for blast_result in blast_results:
                try:
                    pdb_file = 'data/pdb/pdb{}.ent'.format(
                        blast_result.pdb_code)
                    if not os.path.isfile(pdb_file):
                        self.pdb_downloader.retrieve(blast_result.pdb_code)

                    # Silence pdb reader warnings
                    warnings.simplefilter("ignore")
                    angles = torsions.backbone_dihedral_angles(
                        blast_result.pdb_code,
                        blast_result.chain,
                        pdb_file
                    )

                    for hit in blast_result.hits:
                        residue, phi, psi = self.get_central_angles(angles, hit)
                        if residue:
                            torsion_angles['residues'] += residue
                            torsion_angles['psi'].append(phi)
                            torsion_angles['phi'].append(psi)

                except Exception as error:
                    logging.error("Could not download " + blast_result.pdb_code)
                    logging.error(error)

            plot.ramachandran(torsion_angles, fragment)
            secondary_structure = predict_secondary_structure(fragment)
            clusters = cluster_torsion_angles(torsion_angles)
            central_angles = clusters[secondary_structure[self.central]]
            fragment_angles.append(central_angles)
        print(aa_sequence, fragment_angles)
        write_pdb(aa_sequence, fragment_angles, self.central, 'test.pdb')


if __name__ == '__main__':
    if len(sys.argv) > 1:
        app = TerminalApp()
        app.run(sys.argv[1])
    else:
        print('Syntax: cref <aminoacid sequence>')
