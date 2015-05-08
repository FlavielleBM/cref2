import os
import time
import logging
import math

import pandas

from cref import sequence
from cref.sequence.alignment import Blast
from cref.structure import torsions  # , plot
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
    maximum_target_sequences=100,
    blast=dict(
        expect_threshold=1000,
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
        self.params = default_params.copy()
        self.params.update(params)

        self.fragment_size = params['fragment_size']
        self.central = math.floor(self.fragment_size / 2)
        self.blast = Blast(db='data/blastdb/pdbseqres')
        self.ss_predictor = SecondaryStructurePredictor('data/ss.db')
        self.failed_pdbs = []

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
                                 blast_results):
        blast_structures = []

        for blast_result in blast_results:
            pdb_code = blast_result.pdb_code
            chain = blast_result.chain
            try:
                if pdb_code not in self.failed_pdbs:
                    pdb_file = 'data/pdb/{}/pdb{}.ent'.format(
                        pdb_code[1:3], pdb_code)
                    if not os.path.isfile(pdb_file):
                        raise IOError('PDB not available for ' + pdb_code)

                    angles = torsions.backbone_torsion_angles(
                        pdb_file
                    )
                    for hsp in blast_result.hsps:
                        structure = self.get_hsp_structure(
                            pdb_code, chain, fragment, ss, hsp, angles)
                        if structure:
                            blast_structures.append(structure)
            except Exception as e:
                self.failed_pdbs.append(pdb_code)
                logger.warn(e)

        return blast_structures

    def get_secondary_structure(self, aa_sequence, output_dir):
        self.reporter('PREDICTING_SECONDARY_STRUCTURE')
        prediction = self.ss_predictor.sspro(aa_sequence)
        ss = prediction.secondary_structure if prediction else None
        ss_fragments = list(sequence.fragment(ss, self.fragment_size))
        logger.info('SS: ' + ss)
        with open(os.path.join(output_dir, 'secondary_structure.txt'), 'w') as \
                sequence_file:
            sequence_file.write(''.join([ss_eight_to_three(x) for x in ss]))
        ss_fragments = [x.replace('C', '-') for x in ss_fragments]
        return ss_fragments

    def get_angles_for_fragment(self, fragment, ss):
        logger.info('Fragment: ' + fragment)
        self.reporter('RUNNING_BLAST')
        blast_results = self.blast.align(fragment)

        self.reporter('RUNNING_TORSIONS')
        blast_structures = self.get_structures_for_blast(
            fragment, ss, blast_results)
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
        # print(blast_structures[:20].to_string(index=False))
        # plot.ramachandran(blast_structures, fragment, self.central)

        self.reporter('CLUSTERING')
        if len(blast_structures) > 100:
            blast_structures = blast_structures[:100]
        return cluster_torsion_angles(blast_structures, ss[self.central])

    def run(self, aa_sequence, output_dir):
        start_time = time.time()
        self.reporter('STARTED')
        # Aminoacids in the beggining have unknown phi and psi
        dihedral_angles = [(None, None)] * (self.central - 1)
        logger.info('Seq: ' + aa_sequence)

        fragments_ss = self.get_secondary_structure(aa_sequence, output_dir)
        fragments = list(sequence.fragment(aa_sequence, self.fragment_size))
        fragments_len = len(fragments)
        logger.info("Number of fragments: {}".format(fragments_len))

        logger.info('-' * 20)
        for i, fragment in enumerate(fragments):
            logger.info(
                'Progress: {} of {} fragments'.format(i + 1, fragments_len))
            ss = fragments_ss[i]
            angles = self.get_angles_for_fragment(fragment, ss)
            dihedral_angles.append(angles)
            logger.info('-' * 20)

        self.reporter('WRITING_PDB')
        # Amino acids in the end have unbound angles
        dihedral_angles += [(None, None)] * (self.central)
        output_file = os.path.join(output_dir, 'predicted_structure.pdb')
        write_pdb(aa_sequence, dihedral_angles, self.central, output_file)

        elapsed_time = time.time() - start_time
        if elapsed_time > 60:
            logger.info('Prediction took {} minutes'.format(elapsed_time / 60))
        else:
            logger.info('Prediction took {} seconds'.format(elapsed_time))
        return os.path.abspath(output_file)
