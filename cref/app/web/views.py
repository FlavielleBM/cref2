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
    sequence = flask.request.get_json(force=True)['sequence']
    resp = predict_structure.delay(sequence)
    return success({'task_id': resp.id})


@app.route('/predict/<sequence>', methods=['GET'])
def predict_test(sequence):
    resp = predict_structure.delay(sequence)
    return success({'task_id': resp.id})


@app.route('/status/<task_id>')
def status(task_id):
    result = predict_structure.AsyncResult(task_id)
    return success({'state': result.state})


@app.route('/result/<task_id>')
def result(task_id):
    result = predict_structure.AsyncResult(task_id)
    if result.ready():
        return success({'pdb_file': result.get()})
    else:
        return failure('Task is pending')
