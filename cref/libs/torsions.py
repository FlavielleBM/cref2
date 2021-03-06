"""
Wrapper around Torsions

Torsions is a simple program to read a PDB file and calculate backbone
torsion angles.

See http://www.bioinf.org.uk/software/torsions/
"""
import subprocess
import logging

from Bio.PDB.Polypeptide import three_to_one

logger = logging.getLogger('CReF')


def dihedral_angles(pdb_filepath):
    """
    The backbone dihedral angles of proteins are called φ (phi, involving
    the backbone atoms C'-N-Cα-C'), ψ (psi, involving the backbone atoms
    N-Cα-C'-N) and ω (omega, involving the backbone atoms Cα-C'-N-Cα).

    :param pdb_filepath: Path to the pdb file used as input to torsions
    """
    output = subprocess.check_output([
        './libs/torsions',
        pdb_filepath,
    ])
    lines = output.decode('utf-8').split('\n')[2:]  # remove heading
    result = [line.split() for line in lines[:-1]]
    result = list(zip(*result))

    indices = []
    residues = ''
    phi = []
    psi = []
    omega = []

    for i in range(len(result[0])):
        try:
            residues += three_to_one(result[1][i])
        except Exception as e:
            logger.debug(
                'Could not get one letter code for ' + result[1][i])
            logger.debug(e)
        else:
            indices.append(result[0][i])
            phi.append(float(result[2][i]))
            psi.append(float(result[3][i]))
            omega.append(float(result[4][i]))

    return dict(
        residues=residues, phi=phi, psi=psi, omega=omega, indices=indices)
