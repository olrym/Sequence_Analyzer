import unittest
from sequences import sequences


class DNATests(unittest.TestCase):
    def setUp(self):
        self.DNA = sequences.DNA('atgGCAGgCTAa', 'Test DNA')

    def test_attributes(self):
        self.assertEqual(self.DNA.sequence, 'atgGCAGgCTAa', 'Failed to obtain the DNA sequence')
        self.assertEqual(self.DNA.identifier, 'Test DNA', 'Failed to obtain the DNA ID')

    def test_properties_computation(self):
        self.assertEqual(self.DNA.length(), 12, 'DNA length calculation failed')
        self.assertEqual(self.DNA.gc_content(), 50, 'DNA GC-content calculation failed')
        self.assertEqual(self.DNA.mw_ss(), 3972.4, 'DNA MW calculation failed')

    def test_sequence_operations(self):
        self.assertEqual(self.DNA.compl_strand(), 'tacCGTCcGATt', 'Generation of the complementary DNA strand failed')
        self.assertEqual(self.DNA.rev_compl_strand(), 'tTAGcCTGCcat', 'Generation of the reverse complementary DNA strand '
                                                                      'failed')
        self.assertEqual(self.DNA.transcription(), 'augGCAGgCUAa', 'DNA transcription failed')
        self.assertEqual(self.DNA.translation(), 'MAG', 'DNA translation failed')

    def tearDown(self):
        del self.DNA


class RNATests(unittest.TestCase):
    def setUp(self):
        self.RNA = sequences.RNA('augGCAGgCUAa', 'Test RNA')

    def test_attributes(self):
        self.assertEqual(self.RNA.sequence, 'augGCAGgCUAa', 'Failed to obtain the RNA sequence')
        self.assertEqual(self.RNA.identifier, 'Test RNA', 'Failed to obtain the RNA ID')

    def test_properties_computation(self):
        self.assertEqual(self.RNA.length(), 12, 'RNA length calculation failed')
        self.assertEqual(self.RNA.gc_content(), 50, 'RNA GC-content calculation failed')
        self.assertEqual(self.RNA.mw(), 4136.4, 'RNA MW calculation failed')

    def test_sequence_operations(self):
        self.assertEqual(self.RNA.rev_transcription(), 'atgGCAGgCTAa', 'RNA reverse ranscription failed')
        self.assertEqual(self.RNA.translation(), 'MAG', 'RNA translation failed')

    def tearDown(self):
        del self.RNA


class PolypeptideTests(unittest.TestCase):
    def setUp(self):
        self.Polypeptide = sequences.Polypeptide('MAPGEDHYTCRTRYSGSK', 'Test protein')

    def test_attributes(self):
        self.assertEqual(self.Polypeptide.sequence, 'MAPGEDHYTCRTRYSGSK', 'Failed to obtain the protein sequence')
        self.assertEqual(self.Polypeptide.identifier, 'Test protein', 'Failed to obtain the protein ID')

    def test_properties_computation(self):
        self.assertEqual(self.Polypeptide.length(), 18, 'Protein length calculation failed')
        self.assertEqual(self.Polypeptide.mw(), 2039, 'Protein molecular weight calculation failed')
        self.assertEqual(self.Polypeptide.ext_coef(), 2980, 'Protein extinction coefficient calculation failed')
        self.assertEqual(self.Polypeptide.isoelectric_point(), 8.14, 'Protein pI calculation failed')

if __name__ == '__main__':
    unittest.main()