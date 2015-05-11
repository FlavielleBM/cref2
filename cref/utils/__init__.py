import sqlite3
import logging


class Database:
    """
    Wrapper around the database
    """
    def __init__(self, filename):
        self.conn = sqlite3.connect(filename)

    def execute(self, query, args=None):
        try:
            cursor = self.conn.cursor()
            if args:
                cursor.execute(query, args)
            else:
                cursor.execute(query)
            self.conn.commit()
        except sqlite3.OperationalError as e:
            logging.info(e)

    def retrieve(self, query):
        cursor = self.conn.cursor()
        cursor.execute(query)
        return cursor.fetchone()

    def close(self):
        self.conn.close()
