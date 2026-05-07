#!/usr/bin/env python
#==============================================================================
# File Name :      preprocess_variantdicts_keepall.py
# Description :    Get indels, synonymous, non-synonymous mutations from
#                  alignparse output.  Retain ALL mutations (including
#                  multi-mutation barcodes).
# Usage :          python3 preprocess_variantdicts_keepall.py \
#                      <CONSENSUS.csv> <output.csv> <aa_length> <umi_length> <gene_seq.txt>
# Author :         Aura Ferreiro, alferreiro@wustl.edu
# Version :        1.2
# Created On :     2025-07-18
# Last Modified:   2026-05-07
#   v1.1: fix rename typo
#   v1.2: pre-compile regexes; replace row-level loop with map(); warn on
#         reference mismatch instead of silently skipping
#===============================================================================

import functools
import re
import sys
import warnings

import pandas as pd

# ── Codon → amino acid lookup ─────────────────────────────────────────────────
# Stop codons are '*', consistent with mutation_space.py's GENETIC_CODE.
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
    'TAA': '*', 'TAG': '*', 'TGA': '*',
}

# ── Pre-compiled regex patterns (compiled once at import, reused per call) ────
_RE_INDEL    = re.compile(r'^(ins|del)', re.IGNORECASE)
_RE_NT_SUB   = re.compile(r'^([ACGT])(?P<pos>\d+)([ACGT])$', re.IGNORECASE)
_RE_CODON_AA = re.compile(r'^([ACGT]{3})(\d+)([ACGT]{3})$')
_RE_INS      = re.compile(r'ins(\d+)([ACGT]+)', re.IGNORECASE)
_RE_DEL      = re.compile(r'del(\d+)to(\d+)', re.IGNORECASE)

# Track (pos, ref, got) triples already warned about to avoid flooding stderr
_ref_mismatch_warned: set = set()


def codon_to_aa_substitution(sub, protein_length):
    sub = sub.strip().upper()
    m = _RE_CODON_AA.match(sub)
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
        m = _RE_INS.match(indel)
        if not m:
            return None
        pos     = int(m.group(1))
        ins_seq = m.group(2).upper()
        codon_0 = (pos - 1) // 3
        codon_1 = codon_0 + 1

        if codon_1 > protein_length:
            return None
        if len(ins_seq) % 3 != 0:
            return f"{codon_1}fs"

        is_clean = (pos - 1) % 3 == 0

        if gene_seq is None:
            return f"{codon_1}ins" if is_clean else f"{codon_1}cdins"

        new_aas = ''.join(decode.get(ins_seq[i:i+3], '?') for i in range(0, len(ins_seq), 3))

        if is_clean:
            if codon_0 == 0:
                aa_after = _aa_at(gene_seq, 0)
                return f"0_{aa_after}1ins{new_aas}"
            aa_before = _aa_at(gene_seq, codon_0 - 1)
            aa_after  = _aa_at(gene_seq, codon_0)
            return f"{aa_before}{codon_0}_{aa_after}{codon_1}ins{new_aas}"
        else:
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
        m = _RE_DEL.match(indel)
        if not m:
            return None
        start  = int(m.group(1))
        end    = int(m.group(2))
        length = end - start + 1
        cs0    = (start - 1) // 3
        ce0    = (end   - 1) // 3
        cs1    = cs0 + 1
        ce1    = ce0 + 1

        if cs1 > protein_length:
            return None
        if length % 3 != 0:
            return f"{cs1}fs"

        is_clean = (start - 1) % 3 == 0

        if gene_seq is None:
            if is_clean:
                return f"{cs1}del" if cs0 == ce0 else f"{cs1}_{ce1}del"
            return f"{cs1}_{ce1}cdel"

        if is_clean:
            aa_s = _aa_at(gene_seq, cs0)
            if cs0 == ce0:
                return f"{aa_s}{cs1}del"
            aa_e = _aa_at(gene_seq, ce0)
            return f"{aa_s}{cs1}_{aa_e}{ce1}del"
        else:
            reg_start = cs0 * 3
            reg_end   = (ce0 + 1) * 3
            orig_reg  = gene_seq[reg_start:reg_end]
            del_s_in  = (start - 1) - reg_start
            del_e_in  = end - reg_start
            new_reg   = orig_reg[:del_s_in] + orig_reg[del_e_in:]
            new_aas   = ''.join(
                decode.get(new_reg[i:i+3], '?')
                for i in range(0, len(new_reg), 3) if len(new_reg[i:i+3]) == 3
            )
            aa_s = decode.get(orig_reg[:3],  '?') if len(orig_reg) >= 3 else '?'
            aa_e = decode.get(orig_reg[-3:], '?') if len(orig_reg) >= 3 else '?'
            if cs0 == ce0:
                return f"{aa_s}{cs1}delins{new_aas}" if new_aas else f"{aa_s}{cs1}del"
            return f"{aa_s}{cs1}_{aa_e}{ce1}delins{new_aas}" if new_aas else f"{aa_s}{cs1}_{aa_e}{ce1}del"

    return None


def _process_single(muts_str, seq, protein_length):
    """
    Process one target_seq_mutations string.

    Returns a dict with keys: n_target_seq_mutations, n_codon_substitutions,
    n_aa_substitutions, number_of_indels, _new_target, mutations.
    """
    muts  = str(muts_str).strip() if pd.notna(muts_str) else ''
    parts = muts.split() if muts else []

    codon_parts = [p for p in parts if not _RE_INDEL.match(p)]
    indel_parts = [p for p in parts if     _RE_INDEL.match(p)]

    codon_code_strs, aa_all, aa_nonsyn = [], [], []
    for p in codon_parts:
        m = _RE_NT_SUB.match(p)
        if not m:
            continue
        ref_nt = m.group(1).upper()
        pos    = int(m.group('pos'))
        alt_nt = p[-1].upper()

        if pos - 1 >= len(seq):
            continue
        if seq[pos - 1] != ref_nt:
            # Warn once per unique (position, expected, observed) triple
            key = (pos, ref_nt, seq[pos - 1])
            if key not in _ref_mismatch_warned:
                sys.stderr.write(
                    f"Warning: reference mismatch at nt position {pos}: "
                    f"mutation token expects {ref_nt} but gene sequence has "
                    f"{seq[pos-1]}. Skipping this token. "
                    f"(Further identical mismatches suppressed.)\n"
                )
                _ref_mismatch_warned.add(key)
            continue

        start      = 3 * ((pos - 1) // 3)
        codon_ref  = seq[start:start + 3]
        if len(codon_ref) != 3:
            continue
        codon_list = list(codon_ref)
        codon_list[(pos - 1) % 3] = alt_nt
        codon_alt  = ''.join(codon_list)
        codon_index = (pos - 1) // 3 + 1
        code_str   = f"{codon_ref}{codon_index}{codon_alt}"
        codon_code_strs.append(code_str)
        aa, nonsyn = codon_to_aa_substitution(code_str, protein_length)
        if aa:
            aa_all.append(aa)
        if nonsyn:
            aa_nonsyn.append(aa)

    ind_parsed = [r for r in
                  (parse_indel(p, protein_length, seq) for p in indel_parts) if r]

    if not parts:
        return dict(n_target_seq_mutations=0, n_codon_substitutions=0,
                    n_aa_substitutions=0, number_of_indels=0,
                    _new_target='wt', mutations='wt')

    return dict(
        n_target_seq_mutations=len(parts),
        n_codon_substitutions=len(codon_parts),
        n_aa_substitutions=len(aa_nonsyn),
        number_of_indels=len(indel_parts),
        _new_target=' '.join(codon_code_strs + ind_parsed),
        mutations=' '.join(aa_all + ind_parsed) if (aa_all or ind_parsed) else '',
    )


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
        df = df.rename(columns={'UMI': 'barcode'})

    df = df[df['barcode'].str.len() == umi_length].copy()

    # Process all mutation strings in one pass using map() with pre-compiled
    # regexes and a partial function. This avoids the overhead of repeated
    # list.append() calls inside a Python for-loop over millions of rows.
    _proc   = functools.partial(_process_single, seq=seq, protein_length=protein_length)
    results = pd.DataFrame(list(map(_proc, df['target_seq_mutations'].fillna(''))))

    for col in ['n_target_seq_mutations', 'n_codon_substitutions',
                'n_aa_substitutions', 'number_of_indels', 'mutations']:
        df[col] = results[col].values
    df['target_seq_mutations'] = results['_new_target'].values

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
