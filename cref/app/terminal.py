#!/usr/bin/env python

import os
import sys
import logging
import warnings
from math import floor

from cref import sequence
from cref.sequence.alignment import Blast
from cref.structure import PDB, torsions, plot
from cref.structure import write_pdb
from cref.structure.clustering import cluster_torsion_angles
from cref.structure.secondary import (SecondaryStructureBD,
                                      predict_secondary_structure)


class TerminalApp:

    def __init__(self):
        self.pdb_downloader = PDB.PDBDownloader('data/pdb')
        self.blast = Blast(db='tests/blastdb/pdbseqres')
        self.fragment_size = 7
        self.central = floor(self.fragment_size / 2)

    def get_central_angles(self, angles, hsp):
        central = self.central
        pos = central - hsp.query_start + 1
        subject_start = angles['residues'].find(hsp.sbjct)

        if pos >= 0 and subject_start >= 0:
            pos = subject_start + pos
            return (
                angles['residues'][pos], angles['psi'][pos], angles['phi'][pos])
        return ('', None, None)

    def run(self, aa_sequence):
        fragment_angles = []
        phi_psi_table = []
        ss_bd = SecondaryStructureBD()

        for fragment in sequence.fragment(aa_sequence, self.fragment_size):
            torsion_angles = []
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

                    for hsp in blast_result.hsps:
                        residue, phi, psi = self.get_central_angles(angles, hsp)
                        if phi and psi:
                            hsp_seq, hsp_ss = ss_bd.retrieve(
                                blast_result.pdb_code,
                                blast_result.chain
                            )
                            start = hsp.sbjct_start - hsp.query_start
                            end = hsp.sbjct_end + self.fragment_size - (
                                hsp.query_end)
                            phi_psi_table.append((
                                blast_result.pdb_code,
                                blast_result.chain,
                                fragment,
                                hsp.sbjct,
                                hsp_seq[start:end],
                                hsp_ss[start:end],
                                hsp.identities,
                                hsp.score,
                                phi,
                                psi
                            ))
                            structure = hsp_ss[start + self.central]
                            identity = self.fragment_size / hsp.identities
                            torsion_angles.append((residue, phi, psi, identity, structure))

                except Exception as error:
                    logging.error("Could not download " + blast_result.pdb_code)
                    logging.error(error)

            torsion_angles.sort(key=lambda x: x[3], reverse=True)
            plot.ramachandran(torsion_angles[:10], fragment)
            secondary_structure = predict_secondary_structure(fragment)
            clusters = cluster_torsion_angles(torsion_angles)
            central_angles = clusters[secondary_structure[self.central]]
            fragment_angles.append(central_angles)

        phi_psi_table.sort(key=lambda x: x[6], reverse=True)
        print("PDB\tChain\tFrag\tSubj\tSeq\tStruct\tIdent\tScore\tPHI\tPSI")
        for item in phi_psi_table:
            print("%s\t%s\t%s\t%s\t%s\t%s\t%d\t% .2f\t% .2f\t% .2f" % item)

        write_pdb(aa_sequence, fragment_angles, self.central, 'test.pdb')


if __name__ == '__main__':
    if len(sys.argv) > 1:
        app = TerminalApp()
        app.run(sys.argv[1])

        # import cProfile
        # cProfile.run('app.run(sys.argv[1])')
    else:
        print('Syntax: cref <aminoacid sequence>')
