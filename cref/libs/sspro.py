import subprocess
import shutil
import tempfile
import os


class Prediction:
    """
    Hold information regarding predictions
    """
    sequence = None
    secondary_structure = None
    solvent_accessibility = None

    def __repr__(self):
        return "{}\n{}\n{}\n".format(
            self.sequence,
            self.secondary_structure,
            self.solvent_accessibility
        )


def _read_prediction(sequence, prediction_dir):
    prediction = Prediction()
    prediction.sequence = sequence

    with open(os.path.join(prediction_dir, 'tmp.pred.ss8')) as ss8:
        for line in ss8:
            prediction.secondary_structure = line[:-1]
    with open(os.path.join(prediction_dir, 'tmp.pred.acc20')) as acc:
        for line in acc:
            prediction.solvent_accessibility = line[:-1]
    return prediction


def predict(sequence):
    prediction_dir = tempfile.mkdtemp(prefix='ss_prediction_')

    with open(os.path.join(prediction_dir, 'tmp.fasta'), 'w') as fasta_file:
        fasta_file.write('>Seq\n')
        fasta_file.write(sequence)

    subprocess.check_call([
        'libs/SCRATCH-1D_1.0/bin/run_SCRATCH-1D_predictors.sh',
        prediction_dir + '/tmp.fasta',
        prediction_dir + '/tmp.pred',
        '4'
    ])
    prediction = _read_prediction(sequence, prediction_dir)

    try:
        shutil.rmtree(prediction_dir)
    except OSError:
        print('Could not remove prediction dir ' + prediction_dir)

    return prediction
