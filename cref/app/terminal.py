#!/usr/bin/env python

import sys
import os
import logging
from math import floor

import pandas
from Bio.PDB import PDBList

from cref import sequence
from cref.sequence.alignment import Blast
from cref.structure import torsions  # , plot
from cref.structure import write_pdb
from cref.structure.clustering import cluster_torsion_angles
from cref.structure.secondary import SecondaryStructurePredictor
from cref.structure.secondary import ss_eight_to_three


class TerminalApp:

    def __init__(self, fragment_size=5):
        self.pdb_downloader = PDBList(pdb='data/pdb')
        self.blast = Blast(db='data/blastdb/pdbseqres')
        self.fragment_size = fragment_size
        self.central = floor(self.fragment_size / 2)
        self.ss_predictor = SecondaryStructurePredictor('data/ss.db')

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
                    central_ss=hsp_ss[start:end][self.central],
                    identity=identity,
                    score=hsp.score,
                    phi=round(phi, 2),
                    psi=round(psi, 2),
                )

    def get_structures_for_blast(self, fragment, ss,
                                 blast_results, failed_pdbs):
        blast_structures = []

        for blast_result in blast_results:
            pdb_code = blast_result.pdb_code
            chain = blast_result.chain
            try:
                if pdb_code not in failed_pdbs:
                    pdb_file = 'data/pdb/{}/pdb{}.ent'.format(
                        pdb_code[1:3], pdb_code)
                    if not os.path.isfile(pdb_file):
                        raise Exception('PDB not available for ' + pdb_code)

                    angles = torsions.backbone_torsion_angles(
                        pdb_file
                    )
                    for hsp in blast_result.hsps:
                        structure = self.get_hsp_structure(
                            pdb_code, chain, fragment, ss, hsp, angles)
                        if structure:
                            blast_structures.append(structure)
            except Exception as e:
                failed_pdbs.append(pdb_code)
                logging.warn(e)

        return blast_structures

    def run(self, aa_sequence, output_dir, reporter):
        reporter('STARTED')

        # Aminoacids in the beggining have unknown phi and psi
        dihedral_angles = [(None, None)] * (self.central - 1)
        failed_pdbs = []

        print('Seq:', aa_sequence)

        reporter('PREDICTING_SECONDARY_STRUCTURE')
        prediction = self.ss_predictor.sspro(aa_sequence)
        ss = prediction.secondary_structure if prediction else None
        ss_fragments = list(sequence.fragment(ss, self.fragment_size))
        print('SS:', ss)
        with open(os.path.join(output_dir, 'secondary_structure.txt'), 'w') as \
                sequence_file:
            sequence_file.write(''.join([ss_eight_to_three(x) for x in ss]))

        fragments = list(sequence.fragment(aa_sequence, self.fragment_size))
        fragments_len = len(fragments)
        print("Number of fragments:", fragments_len)

        for i, fragment in enumerate(fragments):
            print ("Progress: {} of {} fragments".format(i + 1, fragments_len))
            ss = ss_fragments[i]

            reporter('RUNNING_BLAST')
            blast_results = self.blast.align(fragment)

            reporter('RUNNING_TORSIONS')
            blast_structures = self.get_structures_for_blast(
                fragment, ss, blast_results, failed_pdbs)

            blast_structures = pandas.DataFrame(
                blast_structures,
                columns=[
                    'pdb', 'chain', 'fragment', 'subject', 'subject_full',
                    'fragment_ss', 'subject_ss', 'central_ss', 'identity',
                    'score', 'phi', 'psi'
                ]
            )
            blast_structures = blast_structures.sort(
                ['identity', 'score'], ascending=[0,  0])
            print('-' * 100)
            # print(blast_structures[:20].to_string(index=False))
            # plot.ramachandran(blast_structures, fragment, self.central)

            reporter('CLUSTERING')

            # Use only first 100 for clustering
            if len(blast_structures) > 100:
                blast_structures = blast_structures[:100]

            central_angles = cluster_torsion_angles(
                blast_structures, ss[self.central])
            dihedral_angles.append(central_angles)

        # Amino acids in the end have unbound angles
        dihedral_angles += [(None, None)] * (self.central)
        output_file = os.path.join(output_dir, 'predicted_structure.pdb')
        write_pdb(aa_sequence, dihedral_angles, self.central, output_file)
        return os.path.abspath(output_file)


def terminal_reporter(state):
    print(state[0] + state[1:].lower().replace('_', ' '))


def run_cref(aa_sequence, output_dir, fragment_size=5,
             reporter=terminal_reporter):
    pandas.set_option('display.max_columns', 0)
    pandas.set_option('display.max_rows', 5)
    app = TerminalApp(fragment_size)
    return app.run(aa_sequence, output_dir, reporter)

if __name__ == '__main__':
    if len(sys.argv) == 2:
        run_cref(sys.argv[1], 'predictions/')
    elif len(sys.argv) > 3:
        run_cref(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    else:
        print('Syntax: cref <aminoacid sequence>')
