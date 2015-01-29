import unittest
from cref.structure.torsions import backbone_torsion_angles


class TorsionsTestCase(unittest.TestCase):

    def test_backbone_torsion_angles(self):
        residue = ['GLY', '9999.000', '-86.582', '178.106']
        pdb_filename = 'tests/pdbs/pdb1agt.ent'
        angles = backbone_torsion_angles(pdb_filename)
        self.assertIn(residue, angles)
