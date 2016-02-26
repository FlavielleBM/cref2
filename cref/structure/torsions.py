import logging
import pickle
import sqlite3

from cref.utils import Database

logger = logging.getLogger('CReF')


class TorsionAnglesDB(Database):
    """
    Cache torsion angles calculation
    """
    def create(self):
        parent = super(TorsionAnglesDB, self)
        parent.execute(
            """
            CREATE TABLE IF NOT EXISTS pdb_torsions (
                pdb text, residues text, indices blob,
                phi blob, psi blob, omega blob
            )
            """
        )
        parent.execute(
            """
            CREATE INDEX IF NOT EXISTS IdxPDB ON pdb_torsions(pdb)
            """
        )

    def save(self, pdb_code, residues, indices, phi, psi, omega):
        query = "INSERT INTO pdb_torsions VALUES (?, ?, ?, ?, ?, ?)"
        args = (
            pdb_code.upper(),
            residues,
            sqlite3.Binary(pickle.dumps(indices)),
            sqlite3.Binary(pickle.dumps(phi)),
            sqlite3.Binary(pickle.dumps(psi)),
            sqlite3.Binary(pickle.dumps(omega)),
        )
        super(TorsionAnglesDB, self).execute(query, args)

    def retrieve(self, pdb_code):
        result = super(TorsionAnglesDB, self).retrieve(
            """
            SELECT residues, indices, phi,
                psi, omega FROM pdb_torsions WHERE pdb = '{}'
            """.format(pdb_code.upper())
        )
        if result and len(result) == 3:
            return dict(
                residues=result[0],
                indices=pickle.loads(result[1]),
                phi=pickle.loads(result[2]),
                psi=pickle.loads(result[3]),
                omega=pickle.loads(result[4]),
            )
        else:
            return None


class TorsionsCalculator:

    def __init__(self, cache_db='data/torsions.db'):
        self.torsions_db = TorsionAnglesDB(cache_db)

    def get_angles(self, pdb_code):
        angles = self.torsions_db.retrieve(pdb_code)
        if not angles:
            raise KeyError('Torsion angles for {} not found'.format(pdb_code))
        return angles
