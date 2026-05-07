#!/usr/bin/env python
#==============================================================================
# File Name :      preprocess_variantdicts_keepall.py
# Description :    Get indels, synonymous, non-synonymous mutations from
#                  alignparse output.  Retain ALL mutations (including
#                  multi-mutation barcodes).
# Usage :          python3 preprocess_variantdicts_keepall.py \
#                      <CONSENSUS.csv> <output.csv> <aa_length> <umi_length> <gene_seq.txt>
# Author :         Aura Ferreiro, alferreiro@wustl.edu
# Version :        1.1
# Created On :     2025-07-18
# Last Modified:   2026-05-01  (fix rename typo)
#===============================================================================

import pandas as pd
import re
import sys

decode = {
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TGT': 'C', 'TGC': 'C',
    'GAT': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E',
    'TTT': 'F', 'TTC': 'F',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    'CAT': 'H', 'CAC': 'H',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
    'AAA': 'K', 'AAG': 'K',
    'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATG': 'M',
    'AAT': 'N', 'AAC': 'N',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TGG': 'W',
    'TAT': 'Y', 'TAC': 'Y',
    'TAA': '*', 'TAG': '*', 'TGA': '*'
}


def codon_to_aa_substitution(sub, protein_length):
    sub = sub.strip().upper()
    m = re.match(r'^([ACGT]{3})(\d+)([ACGT]{3})$', sub)
    if not m:
        return (None, False)
    ref_codon, pos_str, alt_codon = m.group(1), m.group(2), m.group(3)
    pos = int(pos_str)
    if pos > protein_length:
        return (None, False)
    ref_aa = decode.get(ref_codon, '?')
    alt_aa = decode.get(alt_codon, '?')
    aa_sub = f"{ref_aa}{pos}{alt_aa}"
    return (aa_sub, ref_aa != alt_aa)


def _aa_at(gene_seq, codon_idx_0):
    """Return single-letter AA for 0-indexed codon position, or '?' if out of range."""
    start = codon_idx_0 * 3
    codon = gene_seq[start:start + 3]
    return decode.get(codon, '?') if len(codon) == 3 else '?'


def parse_indel(indel, protein_length, gene_seq=None):
    """
    Parse an alignparse indel token and return HGVS-style protein notation.

    Insertion tokens:  ins{pos}{ACGT...}    (1-indexed nt position)
    Deletion tokens:   del{start}to{end}    (1-indexed nt positions, inclusive)

    Returns one of:
      '{n}fs'                              frameshift (length not divisible by 3)
      '{AA}{n}del'                         clean whole-codon deletion (1 codon)
      '{AA1}{n1}_{AA2}{n2}del'             clean whole-codon deletion (multi-codon)
      '{AA_b}{n_b}_{AA_a}{n_a}ins{AAs}'   clean whole-codon insertion
      '{AA1}{n1}_{AA2}{n2}delins{AAs}'     codon-disrupting deletion (non-boundary)
      '{AA}{n}delins{AAs}'                 codon-disrupting insertion (non-boundary)

    "Clean" means the indel starts exactly at a codon boundary, so it adds or
    removes whole codons without altering any flanking codon.
    "Codon-disrupting" means the indel is in-frame (length % 3 == 0) but starts
    mid-codon, altering the flanking codons as well.
    """
    if not isinstance(indel, str) or not indel:
        return None
    indel = indel.strip()

    # ── INSERTION ─────────────────────────────────────────────────────────────
    if indel.lower().startswith('ins'):
        m = re.match(r'ins(\d+)([ACGT]+)', indel, flags=re.IGNORECASE)
        if not m:
            return None
        pos     = int(m.group(1))
        ins_seq = m.group(2).upper()
        codon_0 = (pos - 1) // 3          # 0-indexed codon containing the site
        codon_1 = codon_0 + 1             # 1-indexed AA position

        if codon_1 > protein_length:
            return None
        if len(ins_seq) % 3 != 0:
            return f"{codon_1}fs"

        is_clean = (pos - 1) % 3 == 0    # insertion falls exactly between codons

        if gene_seq is None:
            return f"{codon_1}ins" if is_clean else f"{codon_1}cdins"

        new_aas = ''.join(decode.get(ins_seq[i:i+3], '?') for i in range(0, len(ins_seq), 3))

        if is_clean:
            # Clean: inserts whole new codon(s) between codon_0-1 and codon_0
            if codon_0 == 0:
                aa_after = _aa_at(gene_seq, 0)
                return f"0_{aa_after}1ins{new_aas}"
            aa_before = _aa_at(gene_seq, codon_0 - 1)
            aa_after  = _aa_at(gene_seq, codon_0)
            return f"{aa_before}{codon_0}_{aa_after}{codon_1}ins{new_aas}"
        else:
            # Codon-disrupting: insertion splits the codon at codon_0
            region_start = codon_0 * 3
            orig_codon   = gene_seq[region_start:region_start + 3]
            split        = (pos - 1) - region_start
            new_region   = orig_codon[:split] + ins_seq + orig_codon[split:]
            result_aas   = ''.join(
                decode.get(new_region[i:i+3], '?')
                for i in range(0, len(new_region), 3) if len(new_region[i:i+3]) == 3
            )
            orig_aa = decode.get(orig_codon, '?') if len(orig_codon) == 3 else '?'
            return f"{orig_aa}{codon_1}delins{result_aas}"

    # ── DELETION ──────────────────────────────────────────────────────────────
    if indel.lower().startswith('del'):
        m = re.match(r'del(\d+)to(\d+)', indel, flags=re.IGNORECASE)
        if not m:
            return None
        start  = int(m.group(1))
        end    = int(m.group(2))
        length = end - start + 1
        cs0    = (start - 1) // 3        # 0-indexed first affected codon
        ce0    = (end   - 1) // 3        # 0-indexed last affected codon
        cs1    = cs0 + 1                 # 1-indexed AA position
        ce1    = ce0 + 1

        if cs1 > protein_length:
            return None
        if length % 3 != 0:
            return f"{cs1}fs"

        is_clean = (start - 1) % 3 == 0  # deletion begins at a codon boundary

        if gene_seq is None:
            if is_clean:
                return f"{cs1}del" if cs0 == ce0 else f"{cs1}_{ce1}del"
            return f"{cs1}_{ce1}cdel"

        if is_clean:
            # Clean whole-codon deletion
            aa_s = _aa_at(gene_seq, cs0)
            if cs0 == ce0:
                return f"{aa_s}{cs1}del"
            aa_e = _aa_at(gene_seq, ce0)
            return f"{aa_s}{cs1}_{aa_e}{ce1}del"
        else:
            # Codon-disrupting deletion: the deletion spans partial codons at its boundaries
            reg_start = cs0 * 3
            reg_end   = (ce0 + 1) * 3
            orig_reg  = gene_seq[reg_start:reg_end]
            del_s_in  = (start - 1) - reg_start   # deletion start within region (0-indexed)
            del_e_in  = end - reg_start             # deletion end within region (exclusive)
            new_reg   = orig_reg[:del_s_in] + orig_reg[del_e_in:]
            new_aas   = ''.join(
                decode.get(new_reg[i:i+3], '?')
                for i in range(0, len(new_reg), 3) if len(new_reg[i:i+3]) == 3
            )
            aa_s = decode.get(orig_reg[:3],  '?') if len(orig_reg) >= 3 else '?'
            aa_e = decode.get(orig_reg[-3:], '?') if len(orig_reg) >= 3 else '?'
            if cs0 == ce0:
                # Entire deletion within one codon span
                return f"{aa_s}{cs1}delins{new_aas}" if new_aas else f"{aa_s}{cs1}del"
            return f"{aa_s}{cs1}_{aa_e}{ce1}delins{new_aas}" if new_aas else f"{aa_s}{cs1}_{aa_e}{ce1}del"

    return None


def process_variants(input_file, output_file, protein_length, umi_length, gene_seq_file):
    with open(gene_seq_file) as fh:
        seq = ''.join(line.strip() for line in fh if not line.startswith('>')).upper()
    if len(seq) < protein_length * 3:
        sys.stderr.write(
            f"Warning: sequence length {len(seq)} shorter than expected {protein_length*3}\n")

    df = pd.read_csv(input_file)
    if df.columns[0].startswith('Unnamed'):
        df = df.drop(df.columns[0], axis=1)

    # Normalise barcode column name
    if 'UMI_sequence' in df.columns and 'barcode' not in df.columns:
        df = df.rename(columns={'UMI_sequence': 'barcode'})
    if 'barcode_sequence' in df.columns and 'barcode' not in df.columns:
        df = df.rename(columns={'barcode_sequence': 'barcode'})
    if 'UMI' in df.columns and 'barcode' not in df.columns:
        df = df.rename(columns={'UMI': 'barcode'})   # was df.rname() — bug fixed

    df = df[df['barcode'].str.len() == umi_length].copy()

    cols = {
        'n_target_seq_mutations': [],
        'n_codon_substitutions': [],
        'n_aa_substitutions': [],
        'number_of_indels': [],
        'mutations': []
    }
    new_target = []

    for muts in df['target_seq_mutations'].fillna('').astype(str):
        muts = muts.strip()
        parts = muts.split() if muts else []
        cols['n_target_seq_mutations'].append(len(parts))

        codon_parts = [p for p in parts if not re.match(r'^(ins|del)', p, flags=re.IGNORECASE)]
        indel_parts = [p for p in parts if re.match(r'^(ins|del)', p, flags=re.IGNORECASE)]
        cols['n_codon_substitutions'].append(len(codon_parts))
        cols['number_of_indels'].append(len(indel_parts))

        codon_code_strs = []
        aa_all = []
        aa_nonsyn = []
        for p in codon_parts:
            m = re.match(r'^([ACGT])(?P<pos>\d+)([ACGT])$', p, flags=re.IGNORECASE)
            if not m:
                continue
            ref_nt = m.group(1).upper()
            pos = int(m.group('pos'))
            alt_nt = p[-1].upper()
            if pos - 1 < len(seq) and seq[pos - 1] != ref_nt:
                continue
            start = 3 * ((pos - 1) // 3)
            codon_ref = seq[start:start + 3]
            if len(codon_ref) != 3:
                continue
            codon_list = list(codon_ref)
            codon_list[(pos - 1) % 3] = alt_nt
            codon_alt = ''.join(codon_list)
            codon_index = (pos - 1) // 3 + 1
            code_str = f"{codon_ref}{codon_index}{codon_alt}"
            codon_code_strs.append(code_str)
            aa, nonsyn = codon_to_aa_substitution(code_str, protein_length)
            if aa:
                aa_all.append(aa)
            if nonsyn:
                aa_nonsyn.append(aa)
        cols['n_aa_substitutions'].append(len(aa_nonsyn))

        ind_parsed = [parse_indel(p, protein_length, seq) for p in indel_parts]
        ind_parsed = [x for x in ind_parsed if x]

        if cols['n_target_seq_mutations'][-1] == 0:
            new_target.append('wt')
            cols['mutations'].append('wt')
        else:
            new_target.append(' '.join(codon_code_strs + ind_parsed))
            cols['mutations'].append(' '.join(aa_all + ind_parsed) if (aa_all or ind_parsed) else '')

    for k, v in cols.items():
        df[k] = v
    df['target_seq_mutations'] = new_target

    out_cols = ['library', 'barcode', 'variant_call_support', 'target_seq_mutations',
                'n_target_seq_mutations', 'n_codon_substitutions',
                'n_aa_substitutions', 'number_of_indels', 'mutations']
    df.to_csv(output_file, columns=out_cols, index=False)
    print(f"Wrote {len(df)} barcodes to {output_file}")


if __name__ == '__main__':
    if len(sys.argv) != 6:
        print("Usage: python3 preprocess_variantdicts_keepall.py "
              "<variants.csv> <output.csv> <aa_length> <umi_length> <gene_seq.txt>")
        sys.exit(1)
    _, variants_file, output_csv, prot_len, umi_len, gene_seq_file = sys.argv
    process_variants(variants_file, output_csv, int(prot_len), int(umi_len), gene_seq_file)
