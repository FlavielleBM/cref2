#!/usr/bin/env python

import sys
import logging
from math import floor

import pandas
from Bio.PDB import PDBList

from cref import sequence
from cref.sequence.alignment import Blast
from cref.structure import torsions, plot
from cref.structure import write_pdb
from cref.structure.clustering import cluster_torsion_angles
from cref.structure.secondary import SecondaryStructurePredictor


class TerminalApp:

    def __init__(self, fragment_size=5):
        self.pdb_downloader = PDBList(pdb='data/pdb')
        self.blast = Blast(db='data/blastdb/pdbseqres')
        self.fragment_size = fragment_size
        self.central = floor(self.fragment_size / 2)
        self.ss_predictor = SecondaryStructurePredictor()

    def get_central_angles(self, angles, hsp):
        central = self.central
        pos = central - hsp.query_start + 1
        subject_start = angles['residues'].find(hsp.sbjct)

        if pos >= 0 and subject_start >= 0:
            pos = subject_start + pos
            return (
                angles['residues'][pos],
                angles['phi'][pos],
                angles['psi'][pos]
            )
        return ('', None, None)

    def get_hsp_structure(self, pdb, chain, fragment, ss, hsp, angles):
        residue, phi, psi = self.get_central_angles(angles, hsp)
        pdb_dssp_result = self.ss_predictor.pdb_dssp(pdb, chain)
        if phi and psi and pdb_dssp_result:
            hsp_seq, hsp_ss = pdb_dssp_result
            start = hsp.sbjct_start - hsp.query_start
            end = hsp.sbjct_end + self.fragment_size - (
                hsp.query_end)
            identity = round(
                100 * (hsp.identities / self.fragment_size))

            if hsp_ss[start:end]:
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
        ignored_pdbs = []
        with open('data/ignored_pdbs.txt', 'r') as ignored_pdbs_file:
            ignored_pdbs = ignored_pdbs_file.read().splitlines()

        for blast_result in blast_results:
            pdb_code = blast_result.pdb_code
            chain = blast_result.chain
            try:
                if pdb_code not in ignored_pdbs:
                    pdb_file = self.pdb_downloader.retrieve_pdb_file(pdb_code)
                    angles = torsions.backbone_torsion_angles(
                        pdb_file
                    )
                    for hsp in blast_result.hsps:
                        structure = self.get_hsp_structure(
                            pdb_code, chain, fragment, ss, hsp, angles)
                        if structure:
                            blast_structures.append(structure)
            except Exception as e:
                logging.warn(e)
                ignored_pdbs.append(pdb_code)

        with open('data/ignored_pdbs.txt', 'w') as ignored_pdbs_file:
            ignored_pdbs = ignored_pdbs_file.write('\n'.join(ignored_pdbs))
        return blast_structures

    def run(self, aa_sequence, output_file):
        # Aminoacids in the beggining have unknown phi and psi
        dihedral_angles = [(None, None)] * (self.central - 1)

        print('Seq:', aa_sequence)

        prediction = self.ss_predictor.porter(aa_sequence)
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
            print('-' * 100)
            print(blast_structures[:20].to_string(index=False))
            plot.ramachandran(blast_structures, fragment, self.central)
            clusters = cluster_torsion_angles(blast_structures)
            central_angles = clusters[ss[self.central]]
            dihedral_angles.append(central_angles)

        # Amino acids in the end have unbound angles
        dihedral_angles += [(None, None)] * (self.central)
        write_pdb(aa_sequence, dihedral_angles, self.central, output_file)


def run_cref(aa_sequence, output_file='output.pdb', fragment_size=5):
    pandas.set_option('display.max_columns', 0)
    pandas.set_option('display.max_rows', 5)
    app = TerminalApp(fragment_size)
    app.run(aa_sequence, output_file)

if __name__ == '__main__':
    if len(sys.argv) == 2:
        run_cref(sys.argv[1])
    elif len(sys.argv) > 3:
        run_cref(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    else:
        print('Syntax: cref <aminoacid sequence>')
