import os
from celery import Celery

from cref.app import BaseApp

app = Celery(
    'tasks',
    backend='db+sqlite:///data/task_results.db',
    broker='amqp://guest@localhost//'
)


@app.task(bind=True)
def predict_structure(self, sequence, params={}):
    output_dir = os.path.join('predictions', self.request.id)
    os.mkdir(output_dir)

    def reporter(state):
        self.update_state(state=state)

    cref_app = BaseApp(params)
    cref_app.reporter = reporter
    return cref_app.run(
        sequence,
        output_dir,
    )
