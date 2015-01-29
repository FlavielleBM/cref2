import subprocess


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
    lines = output.decode('utf-8').split('\n')[2:]  # from 2 to remove heading
    return [line.split()[1:] for line in lines]
