import subprocess
import logging
import pickle
import sqlite3

from Bio.PDB.Polypeptide import three_to_one

from cref.utils import Database

logger = logging.getLogger('CReF')


def backbone_torsion_angles(pdb_filepath):
    """
    Wrapper around Torsions

    Torsions is a simple program to read a PDB file and calculate backbone
    torsion angles.

    The backbone dihedral angles of proteins are called φ (phi, involving
    the backbone atoms C'-N-Cα-C'), ψ (psi, involving the backbone atoms
    N-Cα-C'-N) and ω (omega, involving the backbone atoms Cα-C'-N-Cα).

    See http://www.bioinf.org.uk/software/torsions/

    :param pdb_filepath: Path to the pdb file used as input to torsions
    """
    output = subprocess.check_output([
        './libs/torsions',
        pdb_filepath,
    ])
    lines = output.decode('utf-8').split('\n')[2:]  # remove heading
    result = [line.split()[1:] for line in lines[:-1]]
    result = list(zip(*result))

    residues = ''
    phi = []
    psi = []

    for i in range(len(result[0])):
        try:
            residues += three_to_one(result[0][i])
        except Exception as e:
            logger.debug(
                'Could not get one letter code for ' + result[0][i])
            logger.debug(e)
        else:
            phi.append(float(result[1][i]))
            psi.append(float(result[2][i]))

    return dict(residues=residues, phi=phi, psi=psi)


class TorsionAnglesDB(Database):
    """
    Cache torsion angles calculation
    """
    def create(self):
        parent = super(TorsionAnglesDB, self)
        parent.execute(
            """
            CREATE TABLE IF NOT EXISTS pdb_torsions (
                pdb text, residues text, phi blob, psi blob
            )
            """
        )
        parent.execute(
            """
            CREATE INDEX IF NOT EXISTS IdxPDB ON pdb_torsions(pdb)
            """
        )

    def save(self, pdb_code, residues, phi, psi):
        query = "INSERT INTO pdb_torsions VALUES (?, ?, ?, ?)"
        args = (
            pdb_code.upper(),
            residues,
            sqlite3.Binary(pickle.dumps(phi)),
            sqlite3.Binary(pickle.dumps(psi))
        )
        super(TorsionAnglesDB, self).execute(query, args)

    def retrieve(self, pdb_code):
        result = super(TorsionAnglesDB, self).retrieve(
            """
            SELECT residues, phi, psi FROM pdb_torsions WHERE pdb = '{}'
            """.format(pdb_code.upper())
        )
        if result and len(result) == 3:
            return dict(
                residues=result[0],
                phi=pickle.loads(result[1]),
                psi=pickle.loads(result[2]),
            )
        else:
            return None


class TorsionsCalculator:

    def __init__(self, cache_db='data/torsions.db'):
        self.torsions_db = TorsionAnglesDB(cache_db)

    def get_angles(self, pdb_code, pdb_filepath):
        angles = self.torsions_db.retrieve(pdb_code)
        if not angles:
            angles = backbone_torsion_angles(pdb_filepath)
            if angles:
                self.torsions_db.save(
                    pdb_code,
                    angles['residues'],
                    angles['phi'],
                    angles['psi']
                )
        return angles
