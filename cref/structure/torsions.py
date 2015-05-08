import warnings
import subprocess
import logging
from math import degrees
import Bio.PDB
from Bio.PDB.Polypeptide import three_to_one

logger = logging.getLogger('CReF')


def backbone_dihedral_angles(pdb_code, pdb_chain, pdb_filepath):
    """
    The backbone dihedral angles of proteins are called φ (phi, involving
    the backbone atoms C'-N-Cα-C'), ψ (psi, involving the backbone atoms
    N-Cα-C'-N) and ω (omega, involving the backbone atoms Cα-C'-N-Cα).

    We *only* obtain phi and psi angles.

    :param pdb_code: Identifier for the pdb
    :param pdb_chain: Identifier for the chain inside the model
    :param pdb_filepath: Path to the pdb file used as input

    :return: List residues and torsion angles
        [residues, phi, psi]
    """
    # Silence pdb reader warnings
    warnings.simplefilter("ignore")
    structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filepath)
    residues = ''
    phi_angles = []
    psi_angles = []
    for model in structure:
        chain = model[pdb_chain]
        polypeptides = Bio.PDB.CaPPBuilder().build_peptides(chain)
        for polypeptide in polypeptides:
            phi_psi = polypeptide.get_phi_psi_list()
            # Convert to degrees
            phi_angles += [degrees(phi) if phi else None for phi, _ in phi_psi]
            psi_angles += [degrees(psi) if psi else None for _, psi in phi_psi]
            # Get one letter code
            residues += ''.join(
                [three_to_one(aa.resname) for aa in polypeptide])
    return dict(residues=residues, phi=phi_angles, psi=psi_angles)


def backbone_torsion_angles(pdb_filepath):
    """
    Wrapper around Torsions

    Torsions is a simple program to read a PDB file and calculate backbone
    torsion angles. It calculates phi, psi and omega and can also calculate
    C-alpha pseudo-torsions.

    See http://www.bioinf.org.uk/software/torsions/

    :param pdb_filepath: Path to the pdb file used as input to torsions
    """
    output = subprocess.check_output([
        './cref/structure/torsions',
        pdb_filepath,
    ])
    lines = output.decode('utf-8').split('\n')[2:]  # remove heading
    result = [line.split()[1:] for line in lines[:-1]]
    result = list(zip(*result))

    residues = ''
    phi = []
    psi = []

    for i in range(len(result[0])):
        try:
            residues += three_to_one(result[0][i])
        except Exception as e:
            logger.debug(
                'Could not get one letter code for ' + result[0][i])
            logger.debug(e)
        else:
            phi.append(float(result[1][i]))
            psi.append(float(result[2][i]))

    return dict(residues=residues, phi=phi, psi=psi)
