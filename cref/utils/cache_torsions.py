#!/usr/bin/env python

import glob
from cref.libs import torsions
from cref.structure.torsions import TorsionAnglesDB


def save_torsions_to_db(db_name='data/torsions.db'):
    db.create()
    pdbs = glob.glob('data/pdb/*/*.ent')
    failed_torsions = []

    for i, pdb_filepath in enumerate(pdbs):
        print(i, pdb_filepath[15:19].upper())
        try:
            angles = torsions.dihedral_angles(pdb_filepath)
            db.save(
                pdb_filepath[15:19],
                angles['residues'],
                angles['indices'],
                angles['phi'],
                angles['psi'],
                angles['omega'],
            )
        except Exception as e:
            print(e)
            failed_torsions.append(pdb_filepath)

    with open('torsion_failed_pdbs.txt', 'w') as failed_torsions_file:
        failed_torsions_file.write(''.join(failed_torsions))


if __name__ == '__main__':
    import sys
    db = TorsionAnglesDB(sys.argv[2])
    db.create()
    save_torsions_to_db(sys.argv[1])
    db.close()
