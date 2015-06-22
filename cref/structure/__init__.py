from peptide import PeptideBuilder
import Bio.PDB


def write_pdb(aa_sequence, fragment_angles, gap_length, filepath):
    """
    Generate pdb file with results

    :param aa_sequence: Amino acid sequence
    :param fragment_angles: Backbone torsion angles
    :param gap_length: Length of the gap at the sequence start and end
    :param filepath: Path to the file to save the pdb
    """
    phi, psi = zip(*fragment_angles)
    structure = PeptideBuilder.make_structure(aa_sequence, phi, psi)
    out = Bio.PDB.PDBIO()
    out.set_structure(structure)
    out.save(filepath)


def rmsd(source, target):
    source = Bio.PDB.parser.get_structure('source', source)
    target = Bio.PDB.parser.get_structure('target', target)

    superimposer = Bio.PDB.Superimposer()
    source_atoms = list(source.get_atoms())
    target_atoms = list(target.get_atoms())[:len(source_atoms)]
    superimposer.set_atoms(source_atoms, target_atoms)
    return superimposer.rms
