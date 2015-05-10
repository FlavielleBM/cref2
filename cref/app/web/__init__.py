from flask import Flask
app = Flask(__name__)

from cref.app import BaseApp


class WebApp(BaseApp):
    """
    Web app, a prediction server
    """
