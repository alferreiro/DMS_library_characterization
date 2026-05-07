#!/usr/bin/env python3
"""
Fix bad samples and add new Fastq11 samples to CombinedFastq.

Part 1: Re-concatenate bad samples using corrected SR009992_2 files.
        Source paths from 260407_TET_DMS_badsamples.csv have two types:
          - Fastq8/SR009992/LIBxxxxxx_*.fastq.gz  → resolved by LIB number in SR009992_2 (scratch)
          - Fastq8/SR009922_TopUp/LIBxxxxxx_*.gz  → same filename, resolved in SR009922_TopUp (scratch)
        Output concatenated files overwrite the bad ones in CombinedFastq.

Part 2: Add new Fastq11/SR011046 samples (1 file pair per sample → symlink, no cat needed).

Part 3: Update combined_samplemap.csv in place.

Usage:
    python fix_bad_samples_and_add_fastq11.py [--dry-run]
"""

import argparse
import csv
import os
import re
import subprocess
import sys

# ─── Configuration ──────────────────────────────────────────────────────────────
SCRATCH_BASE  = "/scratch/rmlab/rmlab_shared3/alferreiro/DMS_PacBio/DantasScreen1_Illumina"
OUTPUT_DIR    = os.path.join(SCRATCH_BASE, "CombinedFastq")
MANIFEST_CSV  = os.path.join(OUTPUT_DIR, "combined_samplemap.csv")

SR009992_2_DIR   = os.path.join(SCRATCH_BASE, "Fastq8/SR009992_2")
SR009922_TOP_DIR = os.path.join(SCRATCH_BASE, "Fastq8/SR009922_TopUp")
FASTQ11_DIR      = os.path.join(SCRATCH_BASE, "Fastq11/SR011046")

BAD_SAMPLES_CSV  = os.path.join(OUTPUT_DIR, "260407_TET_DMS_badsamples.csv")
FASTQ11_SAMPLEMAP = os.path.join(FASTQ11_DIR, "Samplemap2.csv")
# ────────────────────────────────────────────────────────────────────────────────


def extract_lib_number(path):
    """Extract 'LIBxxxxxx' from a file path."""
    m = re.search(r"(LIB\d+)_", os.path.basename(path))
    return m.group(1) if m else None


def build_sr009992_2_lookup():
    """
    Return dict: lib_number -> {'R1': abs_path, 'R2': abs_path}
    Built by scanning filenames in SR009992_2 directory.
    """
    lookup = {}
    for fname in os.listdir(SR009992_2_DIR):
        if not fname.endswith(".fastq.gz"):
            continue
        lib = extract_lib_number(fname)
        if lib is None:
            continue
        abs_path = os.path.join(SR009992_2_DIR, fname)
        read = "R1" if "_R1_" in fname else ("R2" if "_R2_" in fname else None)
        if read:
            lookup.setdefault(lib, {})[read] = abs_path
    return lookup


def resolve_source_path(rel_path, sr009992_2_lookup):
    """
    Convert a relative source path from 260407_TET_DMS_badsamples.csv to
    an absolute scratch path using corrected files.
    """
    if "SR009992/" in rel_path:
        lib = extract_lib_number(rel_path)
        if lib not in sr009992_2_lookup:
            raise ValueError(f"LIB {lib} not found in SR009992_2 — cannot resolve: {rel_path}")
        read = "R1" if "_R1_" in rel_path else "R2"
        if read not in sr009992_2_lookup[lib]:
            raise ValueError(f"No {read} file for {lib} in SR009992_2")
        return sr009992_2_lookup[lib][read]
    elif "SR009922_TopUp/" in rel_path:
        fname = os.path.basename(rel_path)
        abs_path = os.path.join(SR009922_TOP_DIR, fname)
        if not os.path.isfile(abs_path):
            raise FileNotFoundError(f"SR009922_TopUp file missing in scratch: {abs_path}")
        return abs_path
    else:
        raise ValueError(f"Unrecognised source path pattern: {rel_path}")


def concatenate(in_files, out_path, dry_run=False):
    missing = [f for f in in_files if not os.path.isfile(f)]
    if missing:
        raise FileNotFoundError(f"Source file(s) not found: {missing}")
    tag = "[DRY RUN] " if dry_run else ""
    print(f"  {tag}cat {len(in_files)} file(s) → {out_path}")
    for f in in_files:
        print(f"    {f}")
    if not dry_run:
        cmd = ["cat"] + in_files
        with open(out_path, "wb") as fout:
            subprocess.run(cmd, stdout=fout, check=True)


def fix_bad_samples(dry_run=False):
    print(f"Loading bad samples from: {BAD_SAMPLES_CSV}")
    sr009992_2_lookup = build_sr009992_2_lookup()
    print(f"SR009992_2 lookup: {len(sr009992_2_lookup)} LIB entries")

    new_rows = {}
    errors = []

    with open(BAD_SAMPLES_CSV, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            sample = row["Sample"]
            r1_rels = [p.strip() for p in row["SourceR1Files"].split(";")]
            r2_rels = [p.strip() for p in row["SourceR2Files"].split(";")]

            print(f"\n  Sample: {sample}")
            try:
                r1_abs = [resolve_source_path(p, sr009992_2_lookup) for p in r1_rels]
                r2_abs = [resolve_source_path(p, sr009992_2_lookup) for p in r2_rels]
            except (ValueError, FileNotFoundError) as e:
                print(f"  ERROR resolving paths: {e}", file=sys.stderr)
                errors.append(str(e))
                continue

            out_r1 = os.path.join(OUTPUT_DIR, f"{sample}_R1.fastq.gz")
            out_r2 = os.path.join(OUTPUT_DIR, f"{sample}_R2.fastq.gz")

            try:
                concatenate(r1_abs, out_r1, dry_run)
                concatenate(r2_abs, out_r2, dry_run)
            except (FileNotFoundError, subprocess.CalledProcessError) as e:
                print(f"  ERROR concatenating: {e}", file=sys.stderr)
                errors.append(str(e))
                continue

            new_rows[sample] = {
                "Library Name": sample,
                "Combined R1": out_r1,
                "Combined R2": out_r2,
                "Num Source Pairs": len(r1_abs),
                "Source R1 Files": ";".join(r1_abs),
                "Source R2 Files": ";".join(r2_abs),
            }

    return new_rows, errors


def add_fastq11_samples(dry_run=False):
    print(f"Loading Fastq11 samples from: {FASTQ11_SAMPLEMAP}")
    new_rows = {}
    errors = []

    with open(FASTQ11_SAMPLEMAP, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            lib_name = row["Library Name"].strip()
            r1_fname = row["FASTQ Path - Read 1"].strip()
            r2_fname = row["FASTQ Path - Read 2"].strip()

            r1_src = os.path.join(FASTQ11_DIR, r1_fname)
            r2_src = os.path.join(FASTQ11_DIR, r2_fname)
            out_r1 = os.path.join(OUTPUT_DIR, f"{lib_name}_R1.fastq.gz")
            out_r2 = os.path.join(OUTPUT_DIR, f"{lib_name}_R2.fastq.gz")

            print(f"\n  Sample: {lib_name}")
            try:
                concatenate([r1_src], out_r1, dry_run)
                concatenate([r2_src], out_r2, dry_run)
            except (FileNotFoundError, subprocess.CalledProcessError) as e:
                msg = str(e)
                print(f"  ERROR: {msg}", file=sys.stderr)
                errors.append(msg)
                continue

            new_rows[lib_name] = {
                "Library Name": lib_name,
                "Combined R1": out_r1,
                "Combined R2": out_r2,
                "Num Source Pairs": 1,
                "Source R1 Files": r1_src,
                "Source R2 Files": r2_src,
            }

    return new_rows, errors


def update_manifest(bad_rows, fastq11_rows, dry_run=False):
    fieldnames = [
        "Library Name", "Combined R1", "Combined R2",
        "Num Source Pairs", "Source R1 Files", "Source R2 Files",
    ]

    existing = {}
    with open(MANIFEST_CSV, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            existing[row["Library Name"]] = {k: row.get(k, "") for k in fieldnames}

    existing.update(bad_rows)
    existing.update(fastq11_rows)

    tag = "[DRY RUN] " if dry_run else ""
    print(f"\n{tag}Writing manifest: {MANIFEST_CSV} ({len(existing)} samples)")

    if not dry_run:
        with open(MANIFEST_CSV, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            for row in sorted(existing.values(), key=lambda r: r["Library Name"]):
                writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true",
                        help="Print what would be done without writing output")
    args = parser.parse_args()

    all_errors = []

    print("=" * 60)
    print("Part 1: Re-processing bad samples (SR009992 → SR009992_2)")
    print("=" * 60)
    bad_rows, errs = fix_bad_samples(dry_run=args.dry_run)
    all_errors.extend(errs)
    print(f"\nProcessed: {len(bad_rows)} bad samples")

    print("\n" + "=" * 60)
    print("Part 2: Adding new Fastq11/SR011046 samples")
    print("=" * 60)
    f11_rows, errs = add_fastq11_samples(dry_run=args.dry_run)
    all_errors.extend(errs)
    print(f"\nProcessed: {len(f11_rows)} new Fastq11 samples")

    print("\n" + "=" * 60)
    print("Part 3: Updating combined_samplemap.csv")
    print("=" * 60)
    update_manifest(bad_rows, f11_rows, dry_run=args.dry_run)

    if all_errors:
        print(f"\n{len(all_errors)} ERROR(s) occurred:", file=sys.stderr)
        for e in all_errors:
            print(f"  {e}", file=sys.stderr)
        sys.exit(1)

    print("\n=== Done ===")


if __name__ == "__main__":
    main()
