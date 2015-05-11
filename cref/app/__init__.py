import os
import time
import logging
import math

import pandas

from cref import sequence
from cref.sequence.alignment import Blast
from cref.structure.torsions import TorsionsCalculator
from cref.structure import write_pdb
from cref.structure.clustering import cluster_torsion_angles
from cref.structure.secondary import SecondaryStructurePredictor
from cref.structure.secondary import ss_eight_to_three

logger = logging.getLogger('CReF')

default_params = dict(
    fragment_size=5,
    number_of_clusters=8,
    exclude=dict(
        identity=0,
        pdbs=[],
    ),
    max_templates=100,
    blast=dict(
        expect_threshold=100000,
        num_alignments=300,
        word_size=2,
        scoring=dict(
            matrix='PAM30',
            gap_costs='ungapped',
        ),
    ),
)


class BaseApp:

    def __init__(self, params):
        self.set_params(params)
        self.blast = Blast(db='data/blastdb/pdbseqres')
        self.ss_predictor = SecondaryStructurePredictor('data/ss.db')
        self.torsions_calculator = TorsionsCalculator('data/torsions.db')
        self.failed_pdbs = []
        self.torsions = {}

    def set_params(self, params):
        self.params = default_params.copy()
        self.params.update(params)
        self.fragment_size = self.params['fragment_size']
        self.central = math.floor(self.fragment_size / 2)
        self.number_of_clusters = self.params['number_of_clusters']
        self.max_templates = self.params['max_templates']
        self.excluded_pdbs = self.params['exclude']['pdbs']

        self.blast_args = {
            'evalue': self.params['blast']['expect_threshold'],
            'word_size': self.params['blast']['word_size'],
            'matrix': self.params['blast']['scoring']['matrix'],
        }

        if 'num_alignments' in self.params['blast']:
            self.blast_args['num_alignments'] = \
                self.params['blast']['num_alignments']

        gap_costs = self.params['blast']['scoring']['gap_costs']
        if gap_costs == 'ungapped':
            self.blast_args['ungapped'] = True
        else:
            self.blast_args['gap_costs'] = gap_costs

    def reporter(self, state):
        logger.info(state[0] + state[1:].lower().replace('_', ' '))

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

    def get_hsp_structure(self, fragment, ss, hsp, angles):
        residue, phi, psi = self.get_central_angles(angles, hsp)
        pdb_dssp_result = self.ss_predictor.pdb_dssp(hsp.pdb_code, hsp.chain)
        if phi and psi and pdb_dssp_result:
            hsp_seq, hsp_ss = pdb_dssp_result
            start = hsp.sbjct_start - hsp.query_start
            end = hsp.sbjct_end + self.fragment_size - (
                hsp.query_end)
            identity = round(
                100 * (hsp.identities / self.fragment_size))

            if hsp_ss[start:end]:
                return dict(
                    pdb=hsp.pdb_code,
                    chain=hsp.chain,
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

    def get_torsion_angles(self, pdb_code):
        pdb_file = 'data/pdb/{}/pdb{}.ent'.format(pdb_code[1:3], pdb_code)
        if not os.path.isfile(pdb_file):
            raise IOError('PDB not available for ' + pdb_code)

        if pdb_code not in self.torsions:
            angles = self.torsions_calculator.get_angles(pdb_code, pdb_file)
            self.torsions[pdb_code] = angles
        else:
            angles = self.torsions[pdb_code]
        return angles

    def get_structures_for_blast(self, fragment, ss, hsps):
        blast_structures = []
        hsps.sort(key=lambda hsp: (hsp.identities, hsp.score), reverse=True)
        for hsp in hsps:
            if len(blast_structures) < self.max_templates:
                try:
                    pdb_code = hsp.pdb_code
                    if (pdb_code not in self.excluded_pdbs) and (
                            pdb_code not in self.failed_pdbs):
                        angles = self.get_torsion_angles(pdb_code)
                        structure = self.get_hsp_structure(
                            fragment, ss, hsp, angles)
                        if structure:
                            blast_structures.append(structure)
                except Exception as e:
                    self.failed_pdbs.append(pdb_code)
                    logger.warn(e)
            else:
                break

        return blast_structures

    def get_secondary_structure(self, aa_sequence, output_dir):
        self.reporter('PREDICTING_SECONDARY_STRUCTURE')
        prediction = self.ss_predictor.sspro(aa_sequence)
        ss = prediction.secondary_structure if prediction else None
        ss_fragments = list(sequence.fragment(ss, self.fragment_size))
        logger.info('Seq: ' + aa_sequence)
        logger.info('Str: ' + ss)
        with open(os.path.join(output_dir, 'secondary_structure.txt'), 'w') as \
                sequence_file:
            sequence_file.write(''.join([ss_eight_to_three(x) for x in ss]))
        ss_fragments = [x.replace('C', '-') for x in ss_fragments]
        return ss_fragments

    def get_angles_for_fragment(self, fragment, ss):
        logger.info('Fragment: ' + fragment)
        self.reporter('RUNNING_BLAST')
        hsps = self.blast.align(fragment, self.blast_args)

        self.reporter('RUNNING_TORSIONS')
        blast_structures = self.get_structures_for_blast(
            fragment, ss, hsps)
        blast_structures = pandas.DataFrame(
            blast_structures,
            columns=[
                'pdb', 'chain', 'fragment', 'subject', 'subject_full',
                'fragment_ss', 'subject_ss', 'central_ss', 'identity',
                'score', 'phi', 'psi'
            ]
        )
        # print(blast_structures[:10].to_string(index=False))
        # plot.ramachandran(blast_structures, fragment, self.central)

        self.reporter('CLUSTERING')
        if len(blast_structures) > self.max_templates:
            blast_structures = blast_structures[:self.max_templates]
        logger.info('Clustering {} fragments'.format(len(blast_structures)))
        return cluster_torsion_angles(
            blast_structures,
            ss[self.central],
            self.number_of_clusters
        )

    def display_elapsed_time(self, start_time):
        elapsed_time = time.time() - start_time
        if elapsed_time > 60:
            logger.info('Prediction took {} minutes'.format(elapsed_time / 60))
        else:
            logger.info('Prediction took {} seconds'.format(elapsed_time))

    def run(self, aa_sequence, output_dir):
        self.reporter('STARTED')
        start_time = time.time()

        # Aminoacids in the beggining have unknown phi and psi
        dihedral_angles = [(None, None)] * (self.central - 1)

        fragments_ss = self.get_secondary_structure(aa_sequence, output_dir)
        fragments = list(sequence.fragment(aa_sequence, self.fragment_size))
        fragments_len = len(fragments)
        logger.info("Number of fragments: {}".format(fragments_len))
        logger.info('-' * 30)
        for i, fragment in enumerate(fragments):
            logger.info(
                'Progress: {} of {} fragments'.format(i + 1, fragments_len))
            ss = fragments_ss[i]
            angles = self.get_angles_for_fragment(fragment, ss)
            dihedral_angles.append(angles)
            logger.info('-' * 30)

        # Amino acids in the end have unbound angles
        dihedral_angles += [(None, None)] * (self.central)

        self.reporter('WRITING_PDB')
        output_file = os.path.join(output_dir, 'predicted_structure.pdb')
        write_pdb(aa_sequence, dihedral_angles, self.central, output_file)
        self.display_elapsed_time(start_time)
        return os.path.abspath(output_file)
