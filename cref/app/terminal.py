#!/usr/bin/env python
import sys
import logging
import warnings
from math import floor
from collections import OrderedDict
import matplotlib.pyplot as plt
import pandas
from cref import sequence
from cref.sequence.alignment import Blast
from cref.structure import PDB, torsions


class TerminalApp:

    def __init__(self):
        self.pdb_downloader = PDB.PDBDownloader('data/pdb')
        self.blast = Blast(db='tests/blastdb/pdbseqres')
        self.fragment_size = 5
        self.ramachandran_densities = pandas.read_csv(
            'data/rama500-general.data',
            skiprows=6,
            delimiter=' ',
            names=['phi', 'psi', 'value']
        )

    def plot_ramachandran(self):
        densities = self.ramachandran_densities
        fontsize = 18
        ticks = [-180, -90, 0, 90, 180]
        plt.contourf(
            list(OrderedDict.fromkeys(densities['phi'])),
            list(OrderedDict.fromkeys(densities['psi'])),
            densities['value'].reshape(180, 180).T,
            levels=[0, 0.0005, 0.02, 1],
            colors=['#FFFFFF', '#B3E8FF', '#7FD9FF']
        )
        plt.xlabel('$\phi$', fontsize=fontsize)
        plt.ylabel('$\psi$', fontsize=fontsize)
        plt.xticks(ticks)
        plt.yticks(ticks)
        plt.tick_params(direction="out")
        plt.margins(0.05)
        ax = plt.axes()
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.spines['left'].set_smart_bounds(True)
        ax.spines['bottom'].set_smart_bounds(True)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')

    def plot_angles(self, torsion_angles, fragment):
        self.plot_ramachandran()
        plt.title("Ramachandran plot for " + fragment)
        plt.scatter(
            torsion_angles['phi'], torsion_angles['psi'],  cmap="b", marker='.')
        plt.show()

    def get_central_angles(self, angles, hit):
        central = floor(self.fragment_size / 2)
        pos = central - hit['query_start'] + 1
        subject_start = angles['residues'].find(hit['subject'])

        if pos >= 0 and subject_start >= 0 and \
                hit['subject'][central] == hit['query'][central]:
            pos = subject_start + pos
            return (
                angles['residues'][pos], angles['psi'][pos], angles['phi'][pos])
        return ('', 0, 0)

    def run(self, aa_sequence):
        for fragment in sequence.fragment(aa_sequence, self.fragment_size):
            torsion_angles = dict(residues='', phi=[], psi=[])
            blast_results = self.blast.align(fragment)

            for blast_result in blast_results:
                try:
                    self.pdb_downloader.retrieve(blast_result.pdb_code)
                    warnings.simplefilter("ignore")
                    angles = torsions.backbone_dihedral_angles(
                        blast_result.pdb_code,
                        blast_result.chain,
                        'data/pdb/pdb{}.ent'.format(blast_result.pdb_code)
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

            self.plot_angles(torsion_angles, fragment)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        app = TerminalApp()
        app.run(sys.argv[1])
    else:
        print('Syntax: cref <aminoacid sequence>')
