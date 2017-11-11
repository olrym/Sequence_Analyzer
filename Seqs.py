class DNA(object):
    """A class for DNA sequences. Attributes are the DNA sequence itself, an identifier for the sequence, and the
    sequence features."""
    __slots__ = ['_sequence', '_identifier', 'features']

    def __init__(self, sequence, identifier='DNA', features=None):
        if not set(sequence).issubset('ATGCNWSMKRYBDHVZatgcnwsmkrybdhvz'):
            raise ValueError('Given sequence is not DNA')
        self._sequence = sequence
        self._identifier = identifier
        self.features = features

    @property
    def sequence(self):
        """Property for keeping self._sequence private."""
        return self._sequence

    @sequence.setter
    def sequence(self, value):
        """Set self._sequence to the value string if the string contains only the nucleotide abbreviations (normal or
        degenerate)."""
        if not set(value).issubset('ATGCNWSMKRYBDHVZatgcnwsmkrybdhvz'):
            raise ValueError('Given sequence is not DNA')
        self._sequence = value

    @property
    def identifier(self):
        """Property for keeping self._identifier private."""
        return self._identifier

    @identifier.setter
    def identifier(self, value):
        """Set self._identifier to the value string."""
        self._identifier = value

    def length(self):
        """Calculate the length of self._sequence."""
        return len(self.sequence)

    def gc_content(self):
        """Calculate the GC-content of the DNA sequence,  if sequence is not degenerate."""
        if set(self.sequence).issubset('ATGCatgc'):
            return (self.sequence.upper().count('G') + self.sequence.upper().count('C')) / len(self.sequence) * 100
        else:
            raise ValueError('The given sequence is degenerate')

    def mw_ss(self):
        """Calculate the molecular weight of the single-stranded DNA molecule, if sequence is not degenerate.."""
        if set(self.sequence).issubset('ATGCatgc'):
            return self.sequence.count('a') * 491.2 + self.sequence.count('A') * 491.2 + self.sequence.count(
                'c') * 467.2 + self.sequence.count('C') * 467.2 + self.sequence.count('g') * 507.2 + self.sequence.count(
                'G') * 507.2 + self.sequence.count('t') * 482.2 + self.sequence.count('T') * 482.2
        else:
            raise ValueError('The given sequence is degenerate')

    def mw_ds(self):
        """Calculate the molecular weight of the double-stranded DNA molecule, if sequence is not degenerate."""
        if set(self.sequence).issubset('ATGCatgc'):
            return self.mw_ss() * 2
        else:
            raise ValueError('The given sequence is degenerate')


    def compl_strand(self):
        """Return the sequence string complementary to self._sequence, if sequence is not degenerate."""
        compl_rules = {'a': 't', 'g': 'c', 't': 'a', 'c': 'g', 'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G'}
        compl_strand = str()
        for nucleotide in self.sequence:
            compl_strand += compl_rules[nucleotide]
        if set(self.sequence).issubset('ATGCatgc'):
            return compl_strand
        else:
            raise ValueError('The given sequence is degenerate')

    def rev_compl_strand(self):
        """Return the sequence string reverse complementary to self._sequence, if sequence is not degenerate."""
        compl_rules = {'a': 't', 'g': 'c', 't': 'a', 'c': 'g', 'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G'}
        rev_compl_strand = str()
        for nucleotide in self.sequence[::-1]:
            rev_compl_strand += compl_rules[nucleotide]
        if set(self.sequence).issubset('ATGCatgc'):
            return rev_compl_strand
        else:
            raise ValueError('The given sequence is degenerate')

    def transcript(self):
        """Return the RNA sequence string corresponding to self._sequence, if sequence is not degenerate."""
        transcr_rules = {'a': 'a', 'g': 'g', 't': 'u', 'c': 'c', 'A': 'A', 'G': 'G', 'T': 'U', 'C': 'C', 'u': 't', 'U':
                         'T'}
        transcr_strand = str()
        for nucleotide in self.sequence:
            transcr_strand += transcr_rules[nucleotide]
        if set(self.sequence).issubset('ATGCatgc'):
            return transcr_strand
        else:
            raise ValueError('The given sequence is degenerate')

    def translation(self):
        """Return the protein sequence string encoded by self._sequence, if sequence is not degenerate."""
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
        if set(self.sequence).issubset('ATGCatgc'):
            return polypeptide
        else:
            raise ValueError('The given sequence is degenerate')


class RNA(object):
    """A class for RNA sequences. Attributes are the DNA sequence itself, an identifier for the sequence, and the
    sequence features"""
    __slots__ = ['_sequence', '_identifier', 'features']

    def __init__(self, sequence, identifier='RNA', features=None):
        if not set(sequence).issubset('AUGCaugc'):
            raise ValueError('Given sequence is not RNA')
        self._sequence = sequence
        self._identifier = identifier
        self.features = features

    @property
    def sequence(self):
        """Property for keeping self._sequence private."""
        return self._sequence

    @sequence.setter
    def sequence(self, value):
        """Set self._sequence to the value string if the string contains only the nucleotide abbreviations."""
        if not set(value).issubset('AUGCaugc'):
            raise ValueError('Given sequence is not RNA')
        self._sequence = value

    @property
    def identifier(self):
        """Property for keeping  self._identifier private."""
        return self._identifier

    @identifier.setter
    def identifier(self, value):
        """Set self._identifier to the value string."""
        self._identifier = value

    def length(self):
        """Calculate the length of self._sequence."""
        return len(self.sequence)

    def gc_content(self):
        """Calculate the GC-content of the RNA sequence."""
        return (self.sequence.count('G') + self.sequence.count('C')) / len(self.sequence) * 100

    def mw(self):
        """Calculate the molecular weight of the RNA molecule."""
        return self.sequence.count('a') * 507.2 + self.sequence.count('A') * 507.2 + self.sequence.count(
            'c') * 483.2 + self.sequence.count('C') * 483.2 + self.sequence.count('g') * 523.2 + self.sequence.count(
            'G') * 523.2 + self.sequence.count('u') * 484.2 + self.sequence.count('U') * 484.2

    def rev_transcript(self):
        """Return the DNA sequence string corresponding to self._sequence."""
        transcr_rules = {'a': 'a', 'g': 'g', 't': 'u', 'c': 'c', 'A': 'A', 'G': 'G', 'T': 'U', 'C': 'C', 'u': 't', 'U':
                         'T'}
        rev_transcr_strand = str()
        for nucleotide in self.sequence:
            rev_transcr_strand += transcr_rules[nucleotide]
        return rev_transcr_strand

    def translation(self):
        """Return the protein sequence string encoded by self._sequence."""
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
    """A class for protein sequences. Attributes are the DNA sequence itself, an identifier for the sequence, and the
    sequence features"""
    __slots__ = ['_sequence', '_identifier', 'features']

    def __init__(self, sequence, identifier='Polypeptide', features=None):
        if not set(sequence).issubset('ACDEFGHIKLMNPQRSTVWY'):
            raise ValueError('Given sequence is not polypeptide')
        self._sequence = sequence
        self._identifier = identifier
        self.features = features

    @property
    def sequence(self):
        """Property for keeping self._sequence private."""
        return self._sequence

    @sequence.setter
    def sequence(self, value):
        """Set self._sequence to the value string if the string contains only the amino acid residues abbreviations."""
        if not set(value).issubset('ACDEFGHIKLMNPQRSTVWY'):
            raise ValueError('Given sequence is not polypeptide')
        self._sequence = value

    @property
    def identifier(self):
        """Property for keeping self._identifier private."""
        return self._identifier

    @identifier.setter
    def identifier(self, value):
        """Set self._identifier to the value string."""

        self._identifier = value

    def length(self):
        """Calculate the length of self._sequence."""
        return len(self.sequence)

    def mw(self):
        """Calculate the molecular weight of the protein molecule."""
        return self.sequence.count('A') * 71 + self.sequence.count('C') * 103 + self.sequence.count('D') * 115 + \
               self.sequence.count('E') * 129 + self.sequence.count('F') * 147 + self.sequence.count('G') * 57 + \
               self.sequence.count('H') * 137 + self.sequence.count('I') * 113 + self.sequence.count('K') * 128 + \
               self.sequence.count('L') * 113 + self.sequence.count('M') * 131 + self.sequence.count('N') * 114 + \
               self.sequence.count('P') * 97 + self.sequence.count('Q') * 128 + self.sequence.count('R') * 156 + \
               self.sequence.count('S') * 87 + self.sequence.count('T') * 101 + self.sequence.count('V') * 99 + \
               self.sequence.count('W') * 186 + self.sequence.count('Y') * 163

    def exst_coef(self):
        """Calculate the exctinction coefficient of the protein molecule at 280 nm wavelength."""
        return self.sequence.count('Y') * 1490 + self.sequence.count('W') * 5500

    def isoelectric_point(self):
        """Calculate the isoelectric point of the protein molecule.
           Algorithm from http://isoelectric.ovh.org/index.html.
           More details in the paper at DOI 10.1186/s13062-016-0159-9."""
        ph = 0.000
        while 1/(1 + 10**(ph - 9.094)) - 1/(1 + 10**(2.869 - ph)) - self.sequence.count('C')/(1 + 10**(7.555 - ph)) - \
              self.sequence.count('D')/(1 + 10**(3.872 - ph)) - self.sequence.count('E')/(1 + 10**(4.412 - ph)) - \
              self.sequence.count('Y')/(1 + 10**(10.85 - ph)) + self.sequence.count('H')/(1 + 10**(ph - 5.637)) + \
              self.sequence.count('K')/(1 + 10**(ph - 9.052)) + self.sequence.count('R')/(1 + 10**(ph - 12.503)) > 0:
                ph += 0.001
        return ph


def read_dna_fasta(filename):
    """Read the FASTA file with the given filename and return a DNA object with the corresponding sequence and
    identifier."""
    with open(filename) as fasta_dna:
        content = fasta_dna.read()
        return DNA(sequence=content.split('\n', 1)[1].replace('\n', ''), identifier=content.split('\n', 1)[0][1:])


def read_rna_fasta(filename):
    """Read the FASTA file with the given filename and return a RNA object with the corresponding sequence and
    identifier."""
    with open(filename) as fasta_rna:
        content = fasta_rna.read()
        return RNA(sequence=content.split('\n', 1)[1].replace('\n', ''), identifier=content.split('\n', 1)[0][1:])


def read_protein_fasta(filename):
    """Read the FASTA file with the given filename and return a Polypeptide object with the corresponding sequence and
    identifier."""
    with open(filename) as fasta_protein:
        content = fasta_protein.read()
        return Polypeptide(sequence=content.split('\n', 1)[1].replace('\n', ''), identifier=content.split('\n', 1)[0]
               [1:])


def read_dna_multi_fasta(filename):
    """Read the Multi-FASTA file with the given filename and return a list of DNA objects with the corresponding
    sequences and identifiers."""
    dnas = list()
    with open(filename) as multi_fasta_dna:
        content = multi_fasta_dna.read()
        sequences = content.split('>')[1:]
        for sequence in sequences:
            dnas.append(DNA(sequence=sequence.split('\n', 1)[1].replace('\n', ''), identifier=sequence.split('\n', 1)[0]))
    return dnas


def read_rna_multi_fasta(filename):
    """Read the Multi-FASTA file with the given filename and return a list of RNA objects with the corresponding
    sequences and identifiers."""
    RNAs = list()
    with open(filename) as multi_fasta_rna:
        content = multi_fasta_rna.read()
        sequences = content.split('>')[1:]
        for sequence in sequences:
            RNAs.append(RNA(sequence=sequence.split('\n', 1)[1].replace('\n', ''), identifier=sequence.split('\n', 1)[0]))
    return RNAs


def read_protein_multi_fasta(filename):
    """Read the Multi-FASTA file with the given filename and return a list of Polypeptide objects with the
    corresponding sequences and identifiers."""
    Proteins = list()
    with open(filename) as multi_fasta_protein:
        content = multi_fasta_protein.read()
        sequences = content.split('>')[1:]
        for sequence in sequences:
            Proteins.append(Polypeptide(sequence=sequence.split('\n', 1)[1].replace('\n', ''), identifier=sequence.split('\n', 1)[0]))
    return Proteins

def write_dna_fasta(filename, DNA):
    """Write DNA object to a FASTA file"""
    with open(filename, 'w') as fasta_dna:
        fasta_dna.write('>' + DNA.identifier + '\n' + DNA.sequence)

def write_rna_fasta(filename, RNA):
    """Write RNA object to a FASTA file"""
    with open(filename, 'w') as fasta_rna:
        fasta_rna.write('>' + RNA.identifier + '\n' + RNA.sequence)

def write_dna_fasta(filename, Protein):
    """Write Polypeptide object to a FASTA file"""
    with open(filename, 'w') as fasta_protein:
        fasta_protein.write('>' + Protein.identifier + '\n' + Protein.sequence)

def read_fastq(filename, cutoff):
    """Read the FASTQ file with the given filename and return a list of DNA objects with the corresponding
    sequences and identifiers for each read that has all the quality values higher than or equal to the cutoff value."""
    a = 0
    quality_scale = '''!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'''
    Seqs = list()
    with open(filename) as fastq_file:
        line_count = 1
        Sequence = DNA('N')
        for line in fastq_file:
            a += 1
            print(a)
            if line_count == 1:
                Sequence.identifier = line.split('@')[1].rstrip()
                line_count = 2
            elif line_count == 2:
                Sequence.sequence = line.rstrip()
                line_count = 3
            elif line_count == 3:
                line_count = 4
            else:
                for character in line:
                    if character == cutoff:
                        Sequence.sequence = 'N'
                    elif character in quality_scale.split(cutoff)[0]:
                        Sequence.sequence = 'N'
                    else:
                        continue
                if Sequence.sequence != 'N':
                    Seqs.append(Sequence)
                line_count = 1
    return Seqs


# dna_test.identifier = 'myDNA'
# print(dna_test.identifier)
# dna_test.sequence = 'ATgGgacGcAttCTaGC'
# print(dna_test.sequence)
# print(dna_test.length())
# print(dna_test.gc_content())
# print(dna_test.mw_ss())
# print(dna_test.mw_ds())
# print(dna_test.compl_strand())
# print(dna_test.rev_compl_strand())
# print(dna_test.transcript())
# print(dna_test.sequence.upper())
# print(dna_test.translation())
# rna_test = RNA('augc')
# rna_test.identifier = 'myRNA'
# print(rna_test.identifier)
# rna_test.sequence = 'AUgGgacGcAuuCUaGC'
# print(rna_test.sequence)
# print(rna_test.length())
# print(rna_test.gc_content())
# print(rna_test.mw())
# print(rna_test.rev_transcript())
# print(rna_test.translation())


