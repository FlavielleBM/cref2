import os
import statistics

from cref.structure import rmsd
from cref.app.terminal import download_pdb, download_fasta, predict_fasta


pdbs = ['1zdd', '1gab']
runs = 5
fragment_sizes = range(5, 13, 2)
number_of_clusters = range(4, 20, 1)

for pdb in pdbs:
    output_dir = 'predictions/evaluation/{}/'.format(pdb)
    try:
        os.mkdir(output_dir)
    except FileExistsError as e:
        print(e)

    for fragment_size in fragment_sizes:
        fasta_file = output_dir + pdb + '.fasta'
        download_fasta(pdb, fasta_file)
        for n in number_of_clusters:
            rmsds = []
            for run in range(runs):
                params = {
                    'pdb': pdb,
                    'fragment_size': fragment_size,
                    'number_of_clusters': n
                }


                prediction_output = output_dir + str(run)
                os.mkdir(prediction_output)
                output_files = predict_fasta(fasta_file, prediction_output, params)
                predicted_structure = output_files[0]
                filepath = os.path.join(
                    os.path.dirname(predicted_structure),
                    'experimental_structure.pdb'
                )
                experimental_structure = download_pdb(pdb, filepath)

                rmsds.append(rmsd(predicted_structure, experimental_structure))
        print(pdb, fragment_size, n, statistics.mean(rmsds), statistics.pstdev(rmsds))
