
def fragment(sequence, size=5):
    """
    Fragment a string sequence using a sliding window given by size

    :param sequence: String containing the sequence
    :param size: Size of the window

    :return: a fragment of the sequence with the given size
    """
    if size > 0 and len(sequence) >= size:
        for i in range(len(sequence) - size + 1):
            yield sequence[i: i + size]
    else:
        raise ValueError((
            "Cannot create fragments with size {} "
            "in a sequence with length {}"
        ).format(size, len(sequence)))
