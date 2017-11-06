class DNA(object):
    __slots__ = ['_sequence', 'identifier', 'features']

    def __init__(self, sequence, identifier='DNA', features=None):
        if not set(sequence).issubset('ATGCatgc'):
            raise ValueError('Given sequence is not DNA')
        self._sequence = sequence
        self.identifier = identifier
        self.features = features

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, value):
        if not set(value).issubset('ATGCatgc'):
            raise ValueError('Given sequence is not DNA')
        self._sequence = value

    def length(self):
        return len(self.sequence)

    def gc_content(self):
        return (self.sequence.count('G') + self.sequence.count('C')) / len(self.sequence) * 100

    def mw_ss(self):
        return self.sequence.count('a') * 491.2 + self.sequence.count('A') * 491.2 + self.sequence.count(
            'c') * 467.2 + self.sequence.count('C') * 467.2 + self.sequence.count('g') * 507.2 + self.sequence.count(
            'G') * 507.2 + self.sequence.count('t') * 482.2 + self.sequence.count('T') * 482.2

    def mw_ds(self):
        return self.mw_ss() * 2

    def compl_strand(self):
        compl_rules = {'a': 't', 'g': 'c', 't': 'a', 'c': 'g', 'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G'}
        compl_strand = str()
        for nucleotide in self.sequence:
            compl_strand += compl_rules[nucleotide]
        return compl_strand

    def rev_compl_strand(self):
        compl_rules = {'a': 't', 'g': 'c', 't': 'a', 'c': 'g', 'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G'}
        rev_compl_strand = str()
        for nucleotide in self.sequence[::-1]:
            rev_compl_strand += compl_rules[nucleotide]
        return rev_compl_strand

    def transcript(self):
        transcr_rules = {'a': 'a', 'g': 'g', 't': 'u', 'c': 'c', 'A': 'A', 'G': 'G', 'T': 'U', 'C': 'C', 'u': 't', 'U':
                         'T'}
        transcr_strand = str()
        for nucleotide in self.sequence:
            transcr_strand += transcr_rules[nucleotide]
        return transcr_strand

    def translation(self):
        translation_rules = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 's', 'TCA': 'S',
                             'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TAA': 'x', 'TAG': 'x', 'TGT': 'C', 'TGC': 'C',
                             'TGA': 'x', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P',
                             'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                             'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
                             'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N',
                             'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V',
                             'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                             'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
                             'GGG': 'G'}
        sequence = self.sequence.upper()
        codons = []
        polypeptide = str()
        while len(sequence) > 3:
            codons.append(sequence[0:3])
            sequence = sequence[3:]
        for codon in codons:
            polypeptide += translation_rules[codon]
        return polypeptide


class RNA(object):
    __slots__ = ['_sequence', 'identifier', 'features']

    def __init__(self, sequence, identifier='RNA', features=None):
        if not set(sequence).issubset('AUGCaugc'):
            raise ValueError('Given sequence is not RNA')
        self._sequence = sequence
        self.identifier = identifier
        self.features = features

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, value):
        if not set(value).issubset('AUGCaugc'):
            raise ValueError('Given sequence is not RNA')
        self._sequence = value

    def length(self):
        return len(self.sequence)

    def gc_content(self):
        return (self.sequence.count('G') + self.sequence.count('C')) / len(self.sequence) * 100

    def mw(self):
        return self.sequence.count('a') * 507.2 + self.sequence.count('A') * 507.2 + self.sequence.count(
            'c') * 483.2 + self.sequence.count('C') * 483.2 + self.sequence.count('g') * 523.2 + self.sequence.count(
            'G') * 523.2 + self.sequence.count('u') * 484.2 + self.sequence.count('U') * 484.2

    def rev_transcript(self):
        transcr_rules = {'a': 'a', 'g': 'g', 't': 'u', 'c': 'c', 'A': 'A', 'G': 'G', 'T': 'U', 'C': 'C', 'u': 't', 'U':
                         'T'}
        rev_transcr_strand = str()
        for nucleotide in self.sequence:
            rev_transcr_strand += transcr_rules[nucleotide]
        return rev_transcr_strand

    def translation(self):
        translation_rules = {'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'UCU': 'S', 'UCC': 's', 'UCA': 'S',
                             'UCG': 'S', 'UAU': 'Y', 'UAC': 'Y', 'UAA': 'x', 'UAG': 'x', 'UGU': 'C', 'UGC': 'C',
                             'UGA': 'x', 'UGG': 'W', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L', 'CCU': 'P',
                             'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                             'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AUU': 'I', 'AUC': 'I', 'AUA': 'I',
                             'AUG': 'M', 'ACU': 'U', 'ACC': 'U', 'ACA': 'U', 'ACG': 'U', 'AAU': 'N', 'AAC': 'N',
                             'AAA': 'K', 'AAG': 'K', 'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GUU': 'V',
                             'GUC': 'V', 'GUA': 'V', 'GUG': 'V', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                             'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G',
                             'GGG': 'G'}
        sequence = self.sequence.upper()
        codons = []
        polypeptide = str()
        while len(sequence) > 3:
            codons.append(sequence[0:3])
            sequence = sequence[3:]
        for codon in codons:
            polypeptide += translation_rules[codon]
        return polypeptide


class Polypeptide(object):
    __slots__ = ['_sequence', 'identifier', 'features']

    def __init__(self, sequence, identifier='Polypeptide', features=None):
        if not set(sequence).issubset('ACDEFGHIKLMNPQRSTVWY'):
            raise ValueError('Given sequence is not polypeptide')
        self._sequence = sequence
        self.identifier = identifier
        self.features = features

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, value):
        if not set(value).issubset('ACDEFGHIKLMNPQRSTVWY'):
            raise ValueError('Given sequence is not polypeptide')
        self._sequence = value

    def length(self):
        return len(self.sequence)

    def mw(self):
        return self.sequence.count('A') * 71 + self.sequence.count('C') * 103 + self.sequence.count('D') * 115 + \
               self.sequence.count('E') * 129 + self.sequence.count('F') * 147 + self.sequence.count('G') * 57 + \
               self.sequence.count('H') * 137 + self.sequence.count('I') * 113 + self.sequence.count('K') * 128 + \
               self.sequence.count('L') * 113 + self.sequence.count('M') * 131 + self.sequence.count('N') * 114 + \
               self.sequence.count('P') * 97 + self.sequence.count('Q') * 128 + self.sequence.count('R') * 156 + \
               self.sequence.count('S') * 87 + self.sequence.count('T') * 101 + self.sequence.count('V') * 99 + \
               self.sequence.count('W') * 186 + self.sequence.count('Y') * 163

    def exst_coef(self):
        return self.sequence.count('Y') * 1490 + self.sequence.count('W') * 5500

    def isoelectric_point(self):
        """ Algorithm after http://isoelectric.ovh.org/index.html
            more details in paper at DOI 10.1186/s13062-016-0159-9 """
        ph = 0.000
        while 1/(1 + 10**(ph - 9.094)) - 1/(1 + 10**(2.869 - ph)) - self.sequence.count('C')/(1 + 10**(7.555 - ph)) - \
              self.sequence.count('D')/(1 + 10**(3.872 - ph)) - self.sequence.count('E')/(1 + 10**(4.412 - ph)) - \
              self.sequence.count('Y')/(1 + 10**(10.85 - ph)) + self.sequence.count('H')/(1 + 10**(ph - 5.637)) + \
              self.sequence.count('K')/(1 + 10**(ph - 9.052)) + self.sequence.count('R')/(1 + 10**(ph - 12.503)) > 0:
                ph += 0.001
        return ph


def dna_fasta(filename):
    with open(filename) as fasta_dna:
        sequence = fasta_dna.read()
        return DNA(sequence=sequence.split('\n', 1)[1].replace('\n', ''), identifier=sequence.split('\n', 1)[0][1:])


def rna_fasta(filename):
    with open(filename) as fasta_rna:
        sequence = fasta_rna.read()
        return RNA(sequence=sequence.split('\n', 1)[1].replace('\n', ''), identifier=sequence.split('\n', 1)[0][1:])


def protein_fasta(filename):
    with open(filename) as fasta_protein:
        sequence = fasta_protein.read()
        return Polypeptide(sequence=sequence.split('\n', 1)[1].replace('\n', ''), identifier=sequence.split('\n', 1)[0]
               [1:])


dna_test = DNA('atgc', identifier='fastaq1')
dna_test.identifier = 'myDNA'
print(dna_test.identifier)
dna_test.sequence = 'ATgGgacGcAttCTaGC'
print(dna_test.sequence)
print(dna_test.length())
print(dna_test.gc_content())
print(dna_test.mw_ss())
print(dna_test.mw_ds())
print(dna_test.compl_strand())
print(dna_test.rev_compl_strand())
print(dna_test.transcript())
print(dna_test.sequence.upper())
print(dna_test.translation())
rna_test = RNA('augc')
rna_test.identifier = 'myRNA'
print(rna_test.identifier)
rna_test.sequence = 'AUgGgacGcAuuCUaGC'
print(rna_test.sequence)
print(rna_test.length())
print(rna_test.gc_content())
print(rna_test.mw())
print(rna_test.rev_transcript())
print(rna_test.translation())
