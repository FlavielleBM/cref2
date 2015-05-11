#!/usr/bin/env python

import glob
from cref.structure.torsions import TorsionAnglesDB, backbone_torsion_angles


def save_torsions_to_db(db_name='data/torsions.db'):
    db.create()
    pdbs = glob.glob('data/pdb/*/*.ent')
    failed_torsions = []

    for i, pdb_code in enumerate(pdbs):
        print(i, pdb_code[15:19].upper())
        try:
            torsions = backbone_torsion_angles(pdb_code)
        except Exception as e:
            print(e)
            failed_torsions.append(pdb_code)
        db.save(
            pdb_code[15:19],
            torsions['residues'],
            torsions['phi'],
            torsions['psi']
        )

    with open('torsion_failed_pdbs.txt', 'w') as failed_torsions_file:
        failed_torsions_file.write(''.join(failed_torsions))


if __name__ == '__main__':
    import sys
    db = TorsionAnglesDB(sys.argv[2])
    db.create()
    save_torsions_to_db(sys.argv[1])
    db.close()
