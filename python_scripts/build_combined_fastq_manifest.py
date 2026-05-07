#!/usr/bin/env python3
"""
Parse all Samplemap2.csv files from DantasScreen1 Illumina runs,
group R1/R2 files by Library Name (potentially split across lanes/TopUp runs),
concatenate all R1s together and all R2s together per sample, and produce
a combined sample mapping CSV.

Usage:
    python build_combined_fastq_manifest.py [--dry-run]
"""

import argparse
import csv
import os
import subprocess
import sys
from collections import defaultdict

# ─── Configuration ──────────────────────────────────────────────────────────────
BASE_LTS = "/lts/rmlab/rmlab_shared3/alferreiro/DMS_PacBio/DantasScreen1_Illumina"
FASTQ_DIRS = ["Fastq7", "Fastq8", "Fastq9", "Fastq10"]
OUTPUT_DIR = "/scratch/rmlab/rmlab_shared3/alferreiro/DMS_PacBio/DantasScreen1_Illumina/CombinedFastq"
MANIFEST_OUT = os.path.join(OUTPUT_DIR, "combined_samplemap.csv")
# ────────────────────────────────────────────────────────────────────────────────


def collect_file_lists(base, fastq_dirs):
    """
    Walk all SR* subdirectories and collect (r1_path, r2_path, library_name)
    tuples from every Samplemap2.csv.

    Returns:
        lib_files: dict  library_name -> list of (r1_abs_path, r2_abs_path)
        seen_paths: set  of absolute paths we have already encountered
                        (to detect duplicates across sample maps)
    """
    lib_files = defaultdict(list)   # library_name -> [(r1, r2), ...]
    seen_pairs = set()              # (r1, r2) tuples already added
    warnings = []

    for fd in fastq_dirs:
        fd_path = os.path.join(base, fd)
        if not os.path.isdir(fd_path):
            warnings.append(f"WARNING: directory not found: {fd_path}")
            continue

        # Sort so that the base SR dir comes before TopUp variants
        sr_dirs = sorted(
            d for d in os.listdir(fd_path)
            if d.startswith("SR") and os.path.isdir(os.path.join(fd_path, d))
        )

        for sr in sr_dirs:
            sr_path = os.path.join(fd_path, sr)
            smap = os.path.join(sr_path, "Samplemap2.csv")
            if not os.path.isfile(smap):
                warnings.append(f"WARNING: no Samplemap2.csv in {sr_path}")
                continue

            with open(smap, newline="") as fh:
                reader = csv.DictReader(fh)
                for row in reader:
                    r1_name = row["FASTQ Path - Read 1"].strip()
                    r2_name = row["FASTQ Path - Read 2"].strip()
                    lib = row["Library Name"].strip()

                    r1_abs = os.path.join(sr_path, r1_name)
                    r2_abs = os.path.join(sr_path, r2_name)

                    # Validate existence
                    for p in [r1_abs, r2_abs]:
                        if not os.path.isfile(p):
                            warnings.append(f"WARNING: file not found: {p}")

                    pair = (r1_abs, r2_abs)
                    if pair in seen_pairs:
                        warnings.append(
                            f"WARNING: duplicate pair skipped for {lib}: {r1_name}"
                        )
                        continue
                    seen_pairs.add(pair)
                    lib_files[lib].append(pair)

    return lib_files, warnings


def safe_filename(name):
    """Convert a Library Name to a safe filename (replace problematic chars)."""
    return name.replace("/", "_").replace(" ", "_")


def build_manifest_and_concatenate(lib_files, output_dir, dry_run=False):
    """
    For each library, cat all R1 files → <lib>_R1.fastq.gz
                          cat all R2 files → <lib>_R2.fastq.gz
    Write combined_samplemap.csv.
    Returns list of (library_name, combined_r1, combined_r2) tuples.
    """
    os.makedirs(output_dir, exist_ok=True)

    manifest_rows = []
    errors = []

    for lib in sorted(lib_files.keys()):
        pairs = lib_files[lib]
        safe_lib = safe_filename(lib)
        out_r1 = os.path.join(output_dir, f"{safe_lib}_R1.fastq.gz")
        out_r2 = os.path.join(output_dir, f"{safe_lib}_R2.fastq.gz")

        r1_files = [p[0] for p in pairs]
        r2_files = [p[1] for p in pairs]

        n = len(pairs)
        print(f"\n{'[DRY RUN] ' if dry_run else ''}Processing: {lib}  ({n} file pair(s))")
        for r1, r2 in pairs:
            print(f"  R1: {r1}")
            print(f"  R2: {r2}")
        print(f"  → {out_r1}")
        print(f"  → {out_r2}")

        if not dry_run:
            for out_f, in_files in [(out_r1, r1_files), (out_r2, r2_files)]:
                cmd = ["cat"] + in_files
                try:
                    with open(out_f, "wb") as fout:
                        result = subprocess.run(cmd, stdout=fout, check=True)
                except subprocess.CalledProcessError as e:
                    msg = f"ERROR concatenating {lib}: {e}"
                    print(msg, file=sys.stderr)
                    errors.append(msg)
                    continue

        manifest_rows.append({
            "Library Name": lib,
            "Combined R1": out_r1,
            "Combined R2": out_r2,
            "Num Source Pairs": n,
            "Source R1 Files": ";".join(r1_files),
            "Source R2 Files": ";".join(r2_files),
        })

    return manifest_rows, errors


def write_manifest(manifest_rows, manifest_path, dry_run=False):
    if dry_run:
        print(f"\n[DRY RUN] Would write manifest to: {manifest_path}")
        return
    fieldnames = [
        "Library Name", "Combined R1", "Combined R2",
        "Num Source Pairs", "Source R1 Files", "Source R2 Files",
    ]
    with open(manifest_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(f"\nManifest written to: {manifest_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true",
                        help="Print what would be done without writing output")
    args = parser.parse_args()

    print("=== Collecting file lists from sample maps ===")
    lib_files, warnings = collect_file_lists(BASE_LTS, FASTQ_DIRS)

    if warnings:
        print("\n--- Warnings ---")
        for w in warnings:
            print(w)

    print(f"\nFound {len(lib_files)} unique Library Names across all sample maps.")
    total_pairs = sum(len(v) for v in lib_files.values())
    multi = {k: v for k, v in lib_files.items() if len(v) > 1}
    print(f"Total source file pairs: {total_pairs}")
    print(f"Libraries with >1 source pair (will be concatenated): {len(multi)}")

    print("\n=== Building combined FASTQs ===")
    manifest_rows, errors = build_manifest_and_concatenate(
        lib_files, OUTPUT_DIR, dry_run=args.dry_run
    )

    write_manifest(manifest_rows, MANIFEST_OUT, dry_run=args.dry_run)

    if errors:
        print(f"\n{len(errors)} ERROR(s) occurred:", file=sys.stderr)
        for e in errors:
            print(e, file=sys.stderr)
        sys.exit(1)

    print("\n=== Done ===")


if __name__ == "__main__":
    main()
