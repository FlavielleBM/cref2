import os
import subprocess
import sys
import tempfile

from peptide import PeptideBuilder
from Bio.PDB import Polypeptide, PDBIO

from cref.libs import torsions


def write_pdb_peptide(output_filepath, sequence, phi, psi, omega):
    structure = PeptideBuilder.make_structure(
        sequence,
        phi,
        psi,
        omega
    )
    out = PDBIO()
    out.set_structure(structure)
    out.save(output_filepath)


def write_pdb_tleap(output_filepath, sequence, phi, psi):
    leapfile = tempfile.NamedTemporaryFile(mode='w', delete=False)
    leapfile.write('source leaprc.ff14SB\n')
    sequence = ' '.join([Polypeptide.one_to_three(aa) for aa in sequence])
    leapfile.write('system  =  sequence {' + sequence + '}\n')
    for i in range(len(phi)):
        leapfile.write((
            'impose system {%s} '
            '{{"N" "CA" "C" "N" %s} '
            '{"C" "N" "CA" "C" %s}}\n') % (
                i + 1,
                phi[i],
                psi[i],
            )
        )
    leapfile.write('savepdb system ' + output_filepath)
    leapfile.write('\nquit\n')

    leapfile.close()
    subprocess.check_call(['tleap', '-f', leapfile.name])
    # os.remove(leapfile.name)


input_filepath = sys.argv[1]
output = torsions.dihedral_angles(input_filepath)
write_pdb_peptide(
    input_filepath[:-3] + 'peptide.pdb',
    output['residues'],
    output['phi'],
    output['psi'],
    output['omega'],
)
write_pdb_tleap(
    input_filepath[:-3] + 'tleap.pdb',
    output['residues'],
    output['phi'],
    output['psi'],
)
