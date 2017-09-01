params = {
    'fragment_size': 5, # Size of the sliding window on the input sequence, must be an odd number
    'number_of_clusters': 8, # Number of clusters to use when clustering templates
    'exclude': {
        'identity': 0, # (0-100) Exclude sequences with the given identity
        'pdbs': [], # Exclude the given pdbs. The pdbs are given as a string list, e.g. ['1gab', '1zdd']
    },
    'max_templates': 100, # Maximum number of templates to retrieve for each sequence window

    # Blast parameters. See https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=BlastHelp#wordsize
    'blast': {
        'expect_threshold': 900000, # Expect value (E) for saving hits (evalue)
        'number_of_alignments': 500, # Show alignments for this number of database sequences (num_alignments)
        'word_size': 2, # Length of initial exact match

        # See: https://www.ncbi.nlm.nih.gov/blast/html/sub_matrix.html
        'scoring': {
            'matrix': 'PAM30',
            'gap_costs': 'ungapped', # If gap costs are provided, they must be a string "gapopen gapextend", e.g. "9 1"
        }
    }
}
