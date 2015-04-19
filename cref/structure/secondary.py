import logging
import sqlite3
import porter_paleale


class Database:
    """
    Wrapper around the database
    """
    def __init__(self, filename='ss.db'):
        self.conn = sqlite3.connect(filename)

    def execute(self, query):
        try:
            cursor = self.conn.cursor()
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


class SecondaryStructureDB(Database):
    """
    Wrapper around the database for secondary structures
    """

    def create(self):
        parent = super(SecondaryStructureDB, self)
        parent.execute(
            """
            CREATE TABLE IF NOT EXISTS pdb_ss (
                pdb text, chain text,
                sequence text, secondary_structure text
            )
            """
        )
        parent.execute(
            """
            CREATE INDEX IF NOT EXISTS IdxPDB ON pdb_ss(pdb, chain)
            """
        )

    def save(self, pdb, chain, sequence, structure):
        super(SecondaryStructureDB, self).execute(
            "INSERT INTO pdb_ss VALUES ('{}', '{}', '{}', '{}')".format(
                pdb, chain, sequence, structure)
        )

    def retrieve(self, pdb, chain):
        return super(SecondaryStructureDB, self).retrieve(
            """
            SELECT sequence, secondary_structure FROM pdb_ss
                WHERE pdb = '{}' AND chain = '{}'
            """.format(pdb.upper(), chain.upper())
        )


class PredictionCache(Database):
    """
    Cache secondary structure predictions
    """

    def create(self):
        parent = super(PredictionCache, self)
        parent.execute(
            """
            CREATE TABLE IF NOT EXISTS predicted_ss (
                sequence text primary key, secondary_structure text,
                solvent_accessibility text
            )
            """
        )
        parent.execute(
            """
            CREATE INDEX IF NOT EXISTS IdxPredicted ON predicted_ss(sequence)
            """
        )

    def save(self, sequence, structure, accessibility):
        super(PredictionCache, self).execute(
            "INSERT INTO predicted_ss VALUES ('{}', '{}', '{}')".format(
                sequence.upper(), structure, accessibility)
        )

    def retrieve(self, sequence):
        return super(PredictionCache, self).retrieve(
            """
            SELECT secondary_structure, solvent_accessibility FROM predicted_ss
                WHERE sequence = '{}'
            """.format(sequence.upper())
        )


# Cache for predicted secondary structures
prediction_cache = PredictionCache('data/ss.db')
prediction_cache.create()


def predict_secondary_structure(sequence):
    """
    Predict the secondary structure of a given sequence

    :param sequence: Amino acid sequence

    :return: Secondary structure and solvent accessibility prediction

        Porter (Secondary Structure):
        H = Helix   (DSSP classes H, G and I)
        E = Strand  (DSSP classes E and B)
        C = Coil    (DSSP classes S, T and .)

        PaleAle (Relative Solvent Accessibility):
        B = very buried      (<=4% accessible)
        b = somewhat buried  (>4% and <=25% accessible)
        e = somewhat exposed (>25% and <=50% accessible)
        E = very exposed     (>50% accessible)


    """
    prediction = prediction_cache.retrieve(sequence)
    if not prediction:
        prediction = porter_paleale.predict(sequence)
        if prediction:
            prediction_cache.save(
                sequence,
                prediction.secondary_structure,
                prediction.solvent_accessibility
            )
    else:
        logging.info('Retrieved cached secondary structure for ' + sequence)
        prediction = porter_paleale.Prediction(
            sequence, prediction[0], prediction[1])
    return prediction
