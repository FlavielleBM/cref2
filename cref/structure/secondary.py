import logging

from cref.libs import sspro
from cref.utils import Database


_ss_eight_to_three = {
    'H': 'H',  # Alpha helix
    'G': 'H',  # 3-10 helix
    'I': 'H',  # Pi helix
    'E': 'E',  # Extended strand
    'B': 'E',  # Beta bridge
    'T': 'C',  # Turn
    'S': 'C',  # Bend
    'C': 'C',  # Other
    '-': 'C',  # Other
}

_ss_similar = {
    'H': ('G', 'I'),
    'G': ('H', 'I'),
    'I': ('H', 'G'),
    'E': ('B',),
    'B': ('E',),
    'T': ('-', 'S'),
    'S': ('-', 'T'),
    'C': ('-', 'T', 'S'),
    '-': ('C', 'T', 'S'),
}


class Prediction:
    """
    Hold information regarding secondary structure predictions
    """
    def __init__(self, sequence, secondary_structure, solvent_accessibility):
        self.sequence = sequence
        self.secondary_structure = secondary_structure
        self.solvent_accessibility = solvent_accessibility

    def __repr__(self):
        return "{}\n{}\n{}\n".format(
            self.sequence,
            self.secondary_structure,
            self.solvent_accessibility
        )


def ss_eight_to_three(structure):
    """
    Convert DSSP secondary structure to three letters code
        H for Helix
        E for Beta strand
        C for Coil
    """
    if structure not in _ss_eight_to_three:
        raise KeyError(structure + ' is not a valid secondary structure')
    return _ss_eight_to_three[structure]


def closest_ss(structure):
    if structure not in _ss_similar:
        raise KeyError(structure + ' is not a valid secondary structure')
    return _ss_similar.get(structure, ('C'))


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


class SecondaryStructureDB(Database):
    """
    Cache secondary structure predictions
    """

    def create(self):
        parent = super(SecondaryStructureDB, self)
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
        super(SecondaryStructureDB, self).execute(
            "INSERT INTO predicted_ss VALUES ('{}', '{}', '{}')".format(
                sequence.upper(), structure, accessibility)
        )

    def retrieve(self, sequence):
        return super(SecondaryStructureDB, self).retrieve(
            """
            SELECT secondary_structure, solvent_accessibility FROM predicted_ss
                WHERE sequence = '{}'
            """.format(sequence.upper())
        )


class SecondaryStructurePredictor:

    def __init__(self, database='data/ss.db'):
        self.prediction_db = SecondaryStructureDB(database)
        self.prediction_db.create()
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

    def _predict(self, sequence, predictor):
        prediction = self.prediction_db.retrieve(sequence)
        if not prediction:
            prediction = predictor(sequence)
            if prediction:
                self.prediction_db.save(
                    sequence,
                    prediction.secondary_structure,
                    prediction.solvent_accessibility
                )
        else:
            logging.info('Read cached secondary structure for ' + sequence)
            prediction = Prediction(sequence, prediction[0], prediction[1])
        return prediction

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
        import porter_paleale
        return self._predict(sequence, porter_paleale.predict)

    def sspro(self, sequence):
        """
        Predict the secondary structure of a given sequence

        :param sequence: Amino acid sequence

        :return Sequence and secondary structure, using the format:
            B = residue in isolated β-bridge
            E = extended strand, participates in β ladder
            G = 3-helix (310 helix)
            I = 5 helix (π-helix)
            T = hydrogen bonded turn
            S = bend
        """
        return self._predict(sequence, sspro.predict)
