import os
from celery import Celery

from cref.app.terminal import run_cref

app = Celery(
    'tasks',
    backend='db+sqlite:///data/task_results.sqlite',
    broker='amqp://guest@localhost//'
)


@app.task(bind=True)
def predict_structure(self, sequence, params={}):
    output_dir = os.path.join('predictions', self.request.id)
    os.mkdir(output_dir)

    def reporter(state):
        self.update_state(state=state)

    return run_cref(
        sequence,
        output_dir,
        reporter=reporter
    )
