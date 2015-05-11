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
    print(flask.request.data)
    params = flask.request.get_json(force=True)
    if 'sequence' not in params:
        return failure('You did not provide an input sequence')
    resp = predict_structure.delay(params['sequence'], params)
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
    return flask.send_from_directory(
        '/home/mchelem/dev/cref2/predictions/',
        filename, as_attachment=True
    )
