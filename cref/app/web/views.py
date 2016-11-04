import os 

import flask
from cref.app.web import app
from cref.app.web.tasks import predict_structure


def success(result):
    return flask.jsonify({
        'status': 'success',
        'retval': result
    })


def failure(reason='Unknown'):
    return flask.jsonify({
        'status': 'failure',
        'reason': reason
    })


@app.route('/predict/', methods=['POST'])
def predict():
    print('\n\n==> Received parameters:', flask.request.data)
    params = flask.request.get_json(force=True)
    if 'sequence' not in params or not params['sequence']:
        print('ERROR: empty input sequence')
        return failure('You did not provide an input sequence')
    sequence = [x for x in params['sequence'] if x in "ACDEFGHIKLMNPQRSTVWYX"]
    resp = predict_structure.delay(''.join(sequence), params)
    return success({'task_id': resp.id})


@app.route('/predict/<sequence>', methods=['GET'])
def predict_test(sequence):
    resp = predict_structure.delay(sequence, {})
    return success({'task_id': resp.id})


@app.route('/status/<task_id>')
def status(task_id):
    result = predict_structure.AsyncResult(task_id)
    return success({'state': result.state})


@app.route('/prediction/<path:filename>')
def download_file(filename):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    dir_path = dir_path[:dir_path.find('cref2') + len('cref2')]
    return flask.send_from_directory(
        os.path.join(dir_path, 'predictions/'),
        filename, as_attachment=True
    )	

