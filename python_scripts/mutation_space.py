#!/usr/bin/env python
#==============================================================================
# File Name :      mutation_space.py
# Description :    Find all possible amino acid mutations for each locus in a
#                  genesequence, assuming only one nucleotide per codon can be
#                  mutated at a time.
# Usage :          python3 mutation_space.py geneseq.txt out.csv
# Author :         Aura Ferreiro, alferreiro@wustl.edu
# Version :        1.1
# Created On :     2024-02-26
# Last Modified:   2026-05-07
#   v1.1: stop codons now encoded as '*' (consistent with preprocess_variantdicts);
#         use .get() to avoid KeyError on ambiguous bases (N, etc.);
#         read_sequence skips FASTA header lines
#===============================================================================

import sys
import pandas as pd

# Stop codons use '*' to match the decode table in preprocess_variantdicts_keepall.py
GENETIC_CODE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

genefile = sys.argv[1]
outfile  = sys.argv[2]


def mutate_codon(codon):
    mutations = []
    nucleotides = ['A', 'T', 'C', 'G']
    for i in range(3):
        for nuc in nucleotides:
            if codon[i] != nuc:
                mutations.append(codon[:i] + nuc + codon[i+1:])
    return mutations


def read_sequence(file_path):
    """Read a nucleotide sequence from a plain-text or FASTA file."""
    with open(file_path, 'r') as fh:
        return ''.join(
            line.strip() for line in fh if not line.startswith('>')
        ).upper()


def analyze_sequence(sequence):
    data = []
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        if len(codon) < 3:
            break
        # Use .get() so ambiguous bases (N, etc.) return '?' rather than KeyError
        original_aa = GENETIC_CODE.get(codon, '?')
        mutations = mutate_codon(codon)
        mutated_aas = {GENETIC_CODE.get(m, '?') for m in mutations} - {'?'}
        data.append({
            'locus': i // 3 + 1,
            'original_aa': original_aa,
            'possible_aa_mutations': ' '.join(sorted(mutated_aas)),
            'num_possible_aa_mutations': len(mutated_aas)
        })
    return pd.DataFrame(data)


def main(file_path, out_path):
    sequence = read_sequence(file_path)
    df = analyze_sequence(sequence)
    df.to_csv(out_path, index=False)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python3 mutation_space.py <geneseq.txt> <out.csv>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
