import unittest
from cref.sequence import fragment

class SequenceTestCase(unittest.TestCase):

    def test_fragment(self):
        fragments = list(fragment('ABCDEFGHIJKLM', 5))
        self.assertEqual(len(fragments), 9)
        self.assertEqual(fragments[0], 'ABCDE')
        self.assertEqual(fragments[-1], 'IJKLM')

    def test_fragment_empty(self):
        fragments = list(fragment(''))
        self.assertEqual(len(fragments), 0)
