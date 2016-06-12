#!/usr/bin/env python

import glob
import os

from cref.libs import torsions
from cref.structure import sidechain
from cref.structure.torsions import TorsionAnglesDB


def save_torsions_to_db(pdbs_filepath, db):
    pdbs = glob.glob(os.path.join(pdbs_filepath, '*/*.ent'))
    print(pdbs)
    failed_torsions = []

    for i, pdb_filepath in enumerate(pdbs):
        print(i, pdb_filepath[-8: -4].upper())
        try:
            save_pdb_torsions_to_db(db, pdb_filepath)
        except Exception as e:
            print(e)
            failed_torsions.append(pdb_filepath)

    with open('torsion_failed_pdbs.txt', 'w') as failed_torsions_file:
        failed_torsions_file.write(''.join(failed_torsions))


def save_pdb_torsions_to_db(db, pdb_filepath):
    angles = torsions.dihedral_angles(pdb_filepath)
    chis = sidechain.chi_angles(pdb_filepath)

    db.save(
        pdb_filepath[-8: -4],
        angles['residues'],
        angles['indices'],
        angles['phi'],
        angles['psi'],
        angles['omega'],
        chis
    )


if __name__ == '__main__':
    import sys
    db_name = sys.argv[2]
    pdbs_filepath = sys.argv[1]

    db = TorsionAnglesDB(db_name)
    db.create()
    save_torsions_to_db(pdbs_filepath, db)
    db.close()
