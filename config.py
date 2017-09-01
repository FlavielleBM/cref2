params = {
    'fragment_size': 5,
    'number_of_clusters': 8,
    'exclude': {
        'identity': 0,
        'pdbs': [],
    },
    'max_templates': 100,
    'blast': {
        'expect_threshold': 900000,
        'number_of_alignments': 500,
        'word_size': 2,
        'scoring': {
            'matrix': 'PAM30',
            'gap_costs': 'ungapped',
        }
    }
}
