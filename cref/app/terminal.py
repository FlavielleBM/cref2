#!/usr/bin/env python

import sys
import logging
from math import floor

import pandas

from cref import sequence
from cref.sequence.alignment import Blast
from cref.structure import PDB, torsions, plot
from cref.structure import write_pdb
from cref.structure.clustering import cluster_torsion_angles
from cref.structure.secondary import (SecondaryStructureDB,
                                      predict_secondary_structure)


class TerminalApp:

    def __init__(self):
        self.pdb_downloader = PDB.PDBDownloader('data/pdb')
        self.blast = Blast(db='tests/blastdb/pdbseqres')
        self.fragment_size = 7
        self.central = floor(self.fragment_size / 2)
        self.ss_db = SecondaryStructureDB()

    def get_central_angles(self, angles, hsp):
        central = self.central
        pos = central - hsp.query_start + 1
        subject_start = angles['residues'].find(hsp.sbjct)

        if pos >= 0 and subject_start >= 0:
            pos = subject_start + pos
            return (
                angles['residues'][pos], angles['phi'][pos], angles['psi'][pos])
        return ('', None, None)

    def get_hsp_structure(self, pdb, chain, fragment, ss, hsp, angles):
        residue, phi, psi = self.get_central_angles(angles, hsp)
        ss_db_result = self.ss_db.retrieve(pdb, chain)
        if phi and psi and ss_db_result:
            hsp_seq, hsp_ss = ss_db_result
            start = hsp.sbjct_start - hsp.query_start
            end = hsp.sbjct_end + self.fragment_size - (
                hsp.query_end)
            identity = round(
                100 * (hsp.identities / self.fragment_size))

            return dict(
                pdb=pdb,
                chain=chain,
                fragment=fragment,
                fragment_ss=ss,
                subject=hsp.sbjct,
                subject_full=hsp_seq[start:end],
                subject_ss=hsp_ss[start:end],
                identity=identity,
                score=hsp.score,
                phi=round(phi, 2),
                psi=round(psi, 2),
            )

    def get_structures_for_blast(self, fragment, ss, blast_results):
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
                        pdb_code, chain, fragment, ss, hsp, angles)
                    if structure:
                        blast_structures.append(structure)

            except Exception as error:
                logging.error("Could not download " + pdb_code)
                logging.error(error)
        return blast_structures

    def run(self, aa_sequence):
        # Aminoacids in the beggining have unknown phi and psi
        dihedral_angles = [(None, None)] * (self.central - 1)

        print('Seq:', aa_sequence)

        prediction = predict_secondary_structure(aa_sequence)
        ss = prediction.secondary_structure if prediction else None
        ss_fragments = list(sequence.fragment(ss, self.fragment_size))
        print('Str:', ss)

        for i, fragment in enumerate(sequence.fragment(
                aa_sequence, self.fragment_size)):
            ss = ss_fragments[i]

            blast_results = self.blast.align(fragment)
            blast_structures = self.get_structures_for_blast(
                fragment, ss, blast_results)

            blast_structures = pandas.DataFrame(
                blast_structures,
                columns=[
                    'pdb', 'chain', 'fragment', 'subject', 'subject_full',
                    'fragment_ss', 'subject_ss', 'identity', 'score',
                    'phi', 'psi'
                ]
            )
            blast_structures = blast_structures.sort(
                ['identity', 'score'], ascending=[0,  0])
            print(blast_structures.to_string(index=False))
            plot.ramachandran(blast_structures, fragment, self.central)
            clusters = cluster_torsion_angles(blast_structures)
            central_angles = clusters[ss[self.central]]
            dihedral_angles.append(central_angles)

        # Amino acids in the end have unbound angles
        dihedral_angles += [(None, None)] * (self.central)
        write_pdb(aa_sequence, dihedral_angles, self.central, 'test.pdb')


def run_cref(args):
    pandas.set_option('display.max_columns', 0)
    pandas.set_option('display.max_rows', 200)
    app = TerminalApp()
    app.run(args)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        run_cref(sys.argv[1])
    else:
        print('Syntax: cref <aminoacid sequence>')
