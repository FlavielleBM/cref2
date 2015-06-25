import os
import statistics

import pandas
import numpy

from cref.structure import rmsd
from cref.app.terminal import download_pdb, run_cref, configure_logger
from cref.app.terminal import download_fasta, read_fasta

pdbs = ['1zdd', '1gb1', '1c5a', '1opd', '1gab']

runs = 10
fragment_sizes = [5, 7, 9, 11, 13, 15]
number_of_clusters = [4, 6, 8, 10, 12, 14, 16, 18]
number_of_templates = [50,  100]

configure_logger()

output_dir = 'predictions/evaluation/'
os.makedirs(output_dir, exist_ok=True)
writer = pandas.ExcelWriter(os.path.join(output_dir, 'results.xlsx'))

for pdb in pdbs:
    results = []
    output_dir = os.path.join('predictions/evaluation/', pdb)
    os.makedirs(output_dir, exist_ok=True)
    fasta_file = os.path.join(output_dir, pdb + '.fasta')
    print('Downloading fasta for', pdb)
    download_fasta(pdb, fasta_file)
    sequence = read_fasta(fasta_file)[0]
    seq = str(sequence.seq).replace('X', '')

    for fragment_size in fragment_sizes:
        for clusters in number_of_clusters:
            for templates in number_of_templates:
                experiment_output = os.path.join(output_dir, '{}_{}_{}'.format(
                    fragment_size, clusters, templates))
                os.makedirs(experiment_output, exist_ok=True)
                rmsds = []
                for run in range(runs):
                    print('>>>', pdb, fragment_size, clusters,  templates, run)
                    params = {
                        'pdb': pdb,
                        'fragment_size': fragment_size,
                        'number_of_clusters': clusters,
                        'max_templates': templates,
                    }
                    prediction_output = os.path.join(
                        experiment_output, str(run))
                    os.makedirs(prediction_output, exist_ok=True)
                    print('Running CReF')
                    predicted_structure = run_cref(
                        seq, prediction_output, params)
                    filepath = os.path.join(
                        os.path.dirname(predicted_structure),
                        'experimental_structure.pdb'
                    )
                    print('Calculating prediction RMSD')
                    experimental_structure = download_pdb(pdb, filepath)
                    rmsds.append(
                        rmsd(predicted_structure, experimental_structure))
                results.append(dict(
                    pdb=pdb,
                    fragment_size=fragment_size,
                    clusters_count=clusters,
                    templates_count=templates,
                    mean_rmsd=statistics.mean(rmsds),
                    sd_rmsd=statistics.pstdev(rmsds),
                    min_rmsd=min(rmsds),
                    min_run=numpy.argmin(rmsds),
                    max_rmsd=max(rmsds),
                    max_run=numpy.argmax(rmsds)
                ))
    results_df = pandas.DataFrame(results)
    results_df.to_excel(os.path.join(output_dir, pdb + '.xlsx'), index=False)
