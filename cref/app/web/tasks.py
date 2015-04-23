from celery import Celery

from cref.app.terminal import run_cref

app = Celery(
    'tasks',
    backend='db+sqlite:///results.sqlite',
    broker='amqp://guest@localhost//'
)


@app.task
def predict_structure(sequence, params={}):
    return run_cref(sequence)
