
def fragment(sequence, size=5):
    """
    Fragment a string sequence using a sliding window given by size

    :param sequence: String containing the sequence
    :param size: Size of the window
    """
    for i in range(len(sequence) - size + 1):
        yield sequence[i: i + size]



