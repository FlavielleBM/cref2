#!/usr/bin/env python

import os
import sys
import logging
import warnings
from math import floor
from pandas import DataFrame

from cref import sequence
from cref.sequence.alignment import Blast
from cref.structure import PDB, torsions, plot
from cref.structure import write_pdb
# from cref.structure.clustering import cluster_torsion_angles
from cref.structure.secondary import SecondaryStructureBD


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
        ss_bd = SecondaryStructureBD()

        for fragment in sequence.fragment(aa_sequence, self.fragment_size):
            cref_results = []
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
                            identity = round(
                                100 * (hsp.identities / self.fragment_size))
                            cref_results.append(dict(
                                pdb=blast_result.pdb_code,
                                chain=blast_result.chain,
                                fragment=fragment,
                                subject=hsp.sbjct,
                                sequence=hsp_seq[start:end],
                                structure=hsp_ss[start:end],
                                identity=identity,
                                score=hsp.score,
                                phi=phi,
                                psi=psi,
                            ))

                except Exception as error:
                    logging.error("Could not download " + blast_result.pdb_code)
                    logging.error(error)

            cref_results = DataFrame(
                cref_results,
                columns=['pdb', 'chain', 'fragment', 'subject', 'sequence',
                         'structure', 'identity', 'score', 'phi', 'psi']
            )
            cref_results = cref_results.sort(
                ['identity', 'score'], ascending=[0,  0])
            print(cref_results)

            plot.ramachandran(cref_results[:30], fragment)
            # secondary_structure = predict_secondary_structure(fragment)
            # clusters = cluster_torsion_angles(torsion_angles)
            # central_angles = clusters[secondary_structure[self.central]]
            # fragment_angles.append(central_angles)

        write_pdb(aa_sequence, fragment_angles, self.central, 'test.pdb')


if __name__ == '__main__':
    if len(sys.argv) > 1:
        app = TerminalApp()
        app.run(sys.argv[1])

        # import cProfile
        # cProfile.run('app.run(sys.argv[1])')
    else:
        print('Syntax: cref <aminoacid sequence>')
