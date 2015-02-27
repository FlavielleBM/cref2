from math import degrees
import Bio.PDB
from Bio.PDB.Polypeptide import three_to_one


def backbone_torsion_angles(pdb_code, pdb_chain, pdb_filepath):
    """
    Wrapper around Torsions

    Torsions is a simple program to read a PDB file and calculate backbone
    torsion angles. It calculates phi, psi and omega and can also calculate
    C-alpha pseudo-torsions.

    See http://www.bioinf.org.uk/software/torsions/

    :param pdb_code: Identifier for the pdb
    :param pdb_chain: Identifier for the chain inside the model
    :param pdb_filepath: Path to the pdb file used as input to torsions

    :return: List of tuples of torsion residues and torsion angles
        (res, phi, psi)
    """
    structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filepath)
    angles = []
    for model in structure:
        chain = model[pdb_chain]
        polypeptides = Bio.PDB.CaPPBuilder().build_peptides(chain)
        for polypeptide in polypeptides:
            phi_psi = polypeptide.get_phi_psi_list()
            # Convert to degrees
            phi = [degrees(phi) if phi else None for phi, _ in phi_psi]
            psi = [degrees(psi) if psi else None for _, psi in phi_psi]
            # Get one letter code
            residues = [three_to_one(aa.get_resname()) for aa in polypeptide]
            angles = list(zip(residues, phi, psi))
    return angles
