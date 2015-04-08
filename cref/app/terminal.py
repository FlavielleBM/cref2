#!/usr/bin/env python

import sys
import logging
from math import floor

import pandas

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
        self.ss_bd = SecondaryStructureBD()

    def get_central_angles(self, angles, hsp):
        central = self.central
        pos = central - hsp.query_start + 1
        subject_start = angles['residues'].find(hsp.sbjct)

        if pos >= 0 and subject_start >= 0:
            pos = subject_start + pos
            return (
                angles['residues'][pos], angles['phi'][pos], angles['psi'][pos])
        return ('', None, None)

    def get_hsp_structure(self, pdb, chain, fragment, hsp, angles):
        residue, phi, psi = self.get_central_angles(angles, hsp)
        if phi and psi:
            hsp_seq, hsp_ss = self.ss_bd.retrieve(pdb, chain)
            start = hsp.sbjct_start - hsp.query_start
            end = hsp.sbjct_end + self.fragment_size - (
                hsp.query_end)
            identity = round(
                100 * (hsp.identities / self.fragment_size))

            return dict(
                pdb=pdb,
                chain=chain,
                fragment=fragment,
                subject=hsp.sbjct,
                sequence=hsp_seq[start:end],
                structure=hsp_ss[start:end],
                identity=identity,
                score=hsp.score,
                phi=phi,
                psi=psi,
            )

    def get_structures_for_blast(self, fragment, blast_results):
        blast_structures = []
        for blast_result in blast_results:
            try:
                pdb_code = blast_result.pdb_code
                chain = blast_result.chain
                pdb_file = self.pdb_downloader.retrieve(pdb_code)

                angles = torsions.backbone_dihedral_angles(
                    pdb_code,
                    chain,
                    pdb_file
                )
                for hsp in blast_result.hsps:
                    structure = self.get_hsp_structure(
                        pdb_code, chain, fragment, hsp, angles)
                    if structure:
                        blast_structures.append(structure)

            except Exception as error:
                logging.error("Could not download " + pdb_code)
                logging.error(error)
        return blast_structures

    def run(self, aa_sequence):
        cref_results = []

        for fragment in sequence.fragment(aa_sequence, self.fragment_size):
            blast_results = self.blast.align(fragment)
            blast_structures = self.get_structures_for_blast(
                fragment, blast_results)

            blast_structures = pandas.DataFrame(
                blast_structures,
                columns=['pdb', 'chain', 'fragment', 'subject', 'sequence',
                         'structure', 'identity', 'score', 'phi', 'psi']
            )
            blast_structures = blast_structures.sort(
                ['identity', 'score'], ascending=[0,  0])
            print(blast_structures)

            plot.ramachandran(blast_structures, fragment, self.central)
            # secondary_structure = predict_secondary_structure(fragment)
            # clusters = cluster_torsion_angles(torsion_angles)
            # central_angles = clusters[secondary_structure[self.central]]
            # fragment_angles.append(central_angles)

        write_pdb(aa_sequence, cref_results, self.central, 'test.pdb')


if __name__ == '__main__':
    if len(sys.argv) > 1:
        pandas.set_option('display.max_columns', 0)
        app = TerminalApp()
        app.run(sys.argv[1])

        # import cProfile
        # cProfile.run('app.run(sys.argv[1])')
    else:
        print('Syntax: cref <aminoacid sequence>')
