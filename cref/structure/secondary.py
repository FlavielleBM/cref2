import sqlite3


class SecondaryStructureBD:
    """
    Wrapper around the database for secondary structures
    """
    def __init__(self, filename='ss.db'):
        self.conn = sqlite3.connect(filename)

    def create(self):
        cursor = self.conn.cursor()
        try:
            cursor.execute(
                '''
                CREATE TABLE pdb_ss (
                    pdb text, chain text,
                    sequence text, secondary_structure text
                )
                '''
            )
            cursor.execute(
                '''
                CREATE INDEX Idx ON pdb_ss(pdb, chain)
                '''
            )
            self.conn.commit()
        except sqlite3.OperationalError as e:
            print(e)

    def save(self, pdb, chain, sequence, structure):
        cursor = self.conn.cursor()
        cursor.execute(
            "INSERT INTO pdb_ss VALUES ('{}', '{}', '{}', '{}')".format(
                pdb, chain, sequence, structure)
        )
        self.conn.commit()

    def close(self):
        self.conn.close()

    def retrieve(self, pdb, chain):
        cursor = self.conn.cursor()
        cursor.execute(
            """
            SELECT sequence, secondary_structure FROM pdb_ss
                WHERE pdb = '{}' AND chain = '{}'
            """.format(pdb.upper(), chain.upper())
        )
        return cursor.fetchone()


def predict_secondary_structure(sequence):
    """
    Predict the secondary structure of a given sequence

    :param sequence: Amino acid sequence

    :return: Secondary structure prediction as a string
        H = helix
        E = Strand
        C =r coil
    """
    # return porter_paleale.predict(sequence)
    return "C" * len(sequence)
