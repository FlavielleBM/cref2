# import porter_paleale


def predict_secondary_structure(sequence):
    """
    Predict the secondary structure of a given sequence

    :param sequence: Amino acid sequence

    :return: Secondary structure prediction as a string
        H = helix
        E = Strand
        C =r coil
    """
    # return porter_paleale.predict(sequence)
    return "C" * len(sequence)


def write_pdb(aa_sequence, fragment_angles, gap_length, filepath):
    """
    Generate pdb file with results

    :param aa_sequence: Amino acid sequence
    :param fragment_angles: Backbone torsion angles
    :param gap_length: Length of the gap at the sequence start and end
    :param filepath: Path to the file to save the pdb
    """
    pass
