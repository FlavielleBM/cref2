import math
import os
import sys

from Bio import PDB

CHI_COUNT = 5
CHI_ATOMS = [
    dict(
        ARG=['N', 'CA', 'CB', 'CG'],
        ASN=['N', 'CA', 'CB', 'CG'],
        ASP=['N', 'CA', 'CB', 'CG'],
        CYS=['N', 'CA', 'CB', 'SG'],
        GLN=['N', 'CA', 'CB', 'CG'],
        GLU=['N', 'CA', 'CB', 'CG'],
        HIS=['N', 'CA', 'CB', 'CG'],
        ILE=['N', 'CA', 'CB', 'CG1'],
        LEU=['N', 'CA', 'CB', 'CG'],
        LYS=['N', 'CA', 'CB', 'CG'],
        MET=['N', 'CA', 'CB', 'CG'],
        PHE=['N', 'CA', 'CB', 'CG'],
        PRO=['N', 'CA', 'CB', 'CG'],
        SER=['N', 'CA', 'CB', 'OG'],
        THR=['N', 'CA', 'CB', 'OG1'],
        TRP=['N', 'CA', 'CB', 'CG'],
        TYR=['N', 'CA', 'CB', 'CG'],
        VAL=['N', 'CA', 'CB', 'CG1'],
    ),
    dict(
        ARG=['CA', 'CB', 'CG', 'CD'],
        ASN=['CA', 'CB', 'CG', 'OD1'],
        ASP=['CA', 'CB', 'CG', 'OD1'],
        GLN=['CA', 'CB', 'CG', 'CD'],
        GLU=['CA', 'CB', 'CG', 'CD'],
        HIS=['CA', 'CB', 'CG', 'ND1'],
        ILE=['CA', 'CB', 'CG1', 'CD1'],
        LEU=['CA', 'CB', 'CG', 'CD1'],
        LYS=['CA', 'CB', 'CG', 'CD'],
        MET=['CA', 'CB', 'CG', 'SD'],
        PHE=['CA', 'CB', 'CG', 'CD1'],
        PRO=['CA', 'CB', 'CG', 'CD'],
        TRP=['CA', 'CB', 'CG', 'CD1'],
        TYR=['CA', 'CB', 'CG', 'CD1'],
    ),
    dict(
        ARG=['CB', 'CG', 'CD', 'NE'],
        GLN=['CB', 'CG', 'CD', 'OE1'],
        GLU=['CB', 'CG', 'CD', 'OE1'],
        LYS=['CB', 'CG', 'CD', 'CE'],
        MET=['CB', 'CG', 'SD', 'CE'],
    ),
    dict(
        ARG=['CG', 'CD', 'NE', 'CZ'],
        LYS=['CG', 'CD', 'CE', 'NZ'],
    ),
    dict(
        ARG=['CD', 'NE', 'CZ', 'NH1'],
    ),
]

def get_sidechain_torsions(filename, model_id):
    parser = PDB.PDBParser()
    structure_id = os.path.splitext(os.path.basename(filename))[0][-4:]
    structure = parser.get_structure(structure_id, filename)
    model = structure.get_list()[model_id]

    torsion_list = []

    for chain in model:
        for res in chain:
            # Skip heteroatoms
            if res.id[0] != " ":
                continue
            res_name = res.resname
            chi_list = [None] * CHI_COUNT
            for i, chi_res in enumerate(CHI_ATOMS):
                if res_name in chi_res:
                    atom_list = chi_res[res_name]
                    vec_atoms = [res[a] for a in atom_list]
                    vectors = [a.get_vector() for a in vec_atoms]
                    angle = PDB.calc_dihedral(*vectors)
                    chi_list[i] = round(math.degrees(angle), 3)

            resi = "{0}{1}".format(res.id[1], res.id[2].strip())
            torsion_list.append([resi, res_name, chi_list])
    return torsion_list

if __name__ == "__main__":
    filename = sys.argv[1]
    model = sys.argv[2] if len(sys.argv) > 2 else 0
    torsion_list = get_sidechain_torsions(filename, model)

    for item in torsion_list:
        print('{:>3}'.format(item[0]), item[1], end=' ')
        for chi in item[2]:
            print('{!s:>8}'.format(chi), end=' ')
        print()
