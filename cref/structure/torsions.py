from math import degrees
import Bio.PDB
from Bio.PDB.Polypeptide import three_to_one


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
