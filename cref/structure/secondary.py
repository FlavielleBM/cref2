import logging
import sqlite3
import porter_paleale


class Database:
    """
    Wrapper around the database
    """
    def __init__(self, filename='data/ss.db'):
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


class PDBSecondaryStructureDB(Database):
    """
    Wrapper around the database for secondary structures
    """

    def create(self):
        parent = super(PDBSecondaryStructureDB, self)
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
        super(PDBSecondaryStructureDB, self).execute(
            "INSERT INTO pdb_ss VALUES ('{}', '{}', '{}', '{}')".format(
                pdb, chain, sequence, structure)
        )

    def retrieve(self, pdb, chain):
        return super(PDBSecondaryStructureDB, self).retrieve(
            """
            SELECT sequence, secondary_structure FROM pdb_ss
                WHERE pdb = '{}' AND chain = '{}'
            """.format(pdb.upper(), chain.upper())
        )


class PorterSecondaryStructureDB(Database):
    """
    Cache secondary structure predictions
    """

    def create(self):
        parent = super(PorterSecondaryStructureDB, self)
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
        super(PorterSecondaryStructureDB, self).execute(
            "INSERT INTO predicted_ss VALUES ('{}', '{}', '{}')".format(
                sequence.upper(), structure, accessibility)
        )

    def retrieve(self, sequence):
        return super(PorterSecondaryStructureDB, self).retrieve(
            """
            SELECT secondary_structure, solvent_accessibility FROM predicted_ss
                WHERE sequence = '{}'
            """.format(sequence.upper())
        )


class SecondaryStructurePredictor:

    def __init__(self, database='data/ss.db'):
        self.porter_db = PorterSecondaryStructureDB(database)
        self.porter_db.create()
        self.pdb_db = PDBSecondaryStructureDB(database)

    def pdb_dssp(self, pdb, chain):
        """
        Get secondary structure obtained from PDB's DSSP

        :param pdb: PDB code
        :param chain: PDB Chain (A, B...)

        :return Sequence and secondary structure, using the format:
            B = residue in isolated β-bridge
            E = extended strand, participates in β ladder
            G = 3-helix (310 helix)
            I = 5 helix (π-helix)
            T = hydrogen bonded turn
            S = bend
        """
        return self.pdb_db.retrieve(pdb, chain)


    def porter(self, sequence):
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
        prediction = self.porter_db.retrieve(sequence)
        if not prediction:
            prediction = porter_paleale.predict(sequence)
            if prediction:
                self.porter_db.save(
                    sequence,
                    prediction.secondary_structure,
                    prediction.solvent_accessibility
                )
        else:
            logging.info('Read cached secondary structure for ' + sequence)
            prediction = porter_paleale.Prediction(
                sequence, prediction[0], prediction[1])
        return prediction
