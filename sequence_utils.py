from collections import namedtuple
import re
import string


complement_translation = string.maketrans('ATCG', 'TAGC')


def reverse_complement(seq):
    return string.translate(seq, complement_translation)[::-1]


def chunks(l, n):
    '''Yield successive n-sized chunks from l.'''
    for i in xrange(0, len(l), n):
        yield l[i:i+n]


def codons(seq):
    return (chunk for chunk in chunks(seq, 3))


def translate(seq):
    return ''.join(CODON_TO_AMINO_ACID[codon] for codon in codons(seq))


def get_transcript_seq(genome_sequence, exons, reverse=False):
    '''Exons are expected to be 1-based closed intervals'''

    exons = sorted(exons, key=lambda e: e.start)
    seq = ''.join(genome_sequence[exon.start - 1:exon.end] for exon in exons)
    if reverse:
        seq = reverse_complement(seq)
    return seq


ORF = namedtuple('ORF', 'start end')

def find_orfs(seq):

    orfs = []

    for match in re.finditer('ATG', seq):
        # Find a start codon
        start = match.start()

        end = 0
        for i, codon in enumerate(codons(seq[start:])):
            if codon in STOP_CODONS:
                # A stop codon was found, set the end position and return the ORF
                end = start + i * 3 + 3
                orfs.append(ORF(start + 1, end))
                break

    return orfs


STOP_CODONS = ('TAA', 'TAG', 'TGA')

CODON_TO_AMINO_ACID = {
    'ATA': 'T', # Isoleucine
    'ATC': 'T', # Isoleucine
    'ATT': 'T', # Isoleucine
    'ATG': 'M', # Methionine
    'ACA': 'T', # Threonine
    'ACC': 'T', # Threonine
    'ACG': 'T', # Threonine
    'ACT': 'T', # Threonine
    'AAC': 'N', # Asparagine
    'AAT': 'N', # Asparagine
    'AAA': 'K', # Lysine
    'AAG': 'K', # Lysine
    'AGC': 'S', # Serine#Valine
    'AGT': 'S', # Serine
    'AGA': 'R', # Arginine
    'AGG': 'R', # Arginine
    'CTA': 'L', # Leucine
    'CTC': 'L', # Leucine
    'CTG': 'L', # Leucine
    'CTT': 'L', # Leucine
    'CCA': 'P', # Proline
    'CAT': 'H', # Histidine
    'CAA': 'Q', # Glutamine
    'CAG': 'Q', # Glutamine
    'CGA': 'R', # Arginine
    'CGC': 'R', # Arginine
    'CGG': 'R', # Arginine
    'CGT': 'R', # Arginine
    'CCC': 'P', # Proline
    'CCG': 'P', # Proline
    'CCT': 'P', # Proline
    'CAC': 'H', # Histidine
    'GTA': 'V', # Valine
    'GTC': 'V', # Valine
    'GTG': 'V', # Valine
    'GTT': 'V', # Valine
    'GCA': 'A', # Alanine
    'GCC': 'A', # Alanine
    'GCG': 'A', # Alanine
    'GCT': 'A', # Alanine
    'GAC': 'D', # Aspartic Acid
    'GAT': 'D', # Aspartic Acid
    'GAA': 'E', # Glutamic Acid
    'GAG': 'E', # Glutamic Acid
    'GGA': 'G', # Glycine
    'GGC': 'G', # Glycine
    'GGG': 'G', # Glycine
    'GGT': 'G', # Glycine
    'TCA': 'S', # Serine
    'TCC': 'S', # Serine
    'TCG': 'S', # Serine
    'TCT': 'S', # Serine
    'TTC': 'F', # Phenylalanine
    'TTT': 'F', # Phenylalanine
    'TTA': 'L', # Leucine
    'TTG': 'L', # Leucine
    'TAC': 'Y', # Tyrosine
    'TAT': 'Y', # Tyrosine
    'TAA': '_', # Stop
    'TAG': '_', # Stop
    'TGC': 'C', # Cysteine
    'TGT': 'C', # Cysteine
    'TGA': '_', # Stop
    'TGG': 'W', # Tryptophan
}
