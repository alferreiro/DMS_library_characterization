#!/usr/bin/env python3
"""
Join bbmerge_perfile.tsv with combined_samplemap.csv on Sample == Library Name.
Output sorted by JoinedPct ascending so poor-merging samples appear first.
"""
import csv

PERFILE  = "/scratch/rmlab/rmlab_shared3/alferreiro/DMS_PacBio/DantasScreen1_Illumina/R1R2merged/bbmerge_perfile.tsv"
SAMPLEMAP = "/scratch/rmlab/rmlab_shared3/alferreiro/DMS_PacBio/DantasScreen1_Illumina/CombinedFastq/combined_samplemap.csv"
OUTFILE  = "/scratch/rmlab/rmlab_shared3/alferreiro/DMS_PacBio/DantasScreen1_Illumina/R1R2merged/bbmerge_with_sources.tsv"

# Load samplemap keyed by Library Name
samplemap = {}
with open(SAMPLEMAP) as f:
    for row in csv.DictReader(f):
        samplemap[row["Library Name"]] = row

# Load perfile and join
rows = []
with open(PERFILE) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        sm = samplemap.get(row["Sample"], {})
        rows.append({
            "Sample":          row["Sample"],
            "Pairs":           int(row["Pairs"]),
            "Joined":          int(row["Joined"]),
            "JoinedPct":       row["JoinedPct"],
            "JoinedPctFloat":  float(row["JoinedPct"].rstrip("%")),
            "NumSourcePairs":  sm.get("Num Source Pairs", "NA"),
            "SourceR1Files":   sm.get("Source R1 Files", "NA"),
            "SourceR2Files":   sm.get("Source R2 Files", "NA"),
        })
        if not sm:
            print(f"WARNING: {row['Sample']} not found in samplemap")

rows.sort(key=lambda r: r["JoinedPctFloat"])

fieldnames = ["Sample", "Pairs", "Joined", "JoinedPct", "NumSourcePairs", "SourceR1Files", "SourceR2Files"]
with open(OUTFILE, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    w.writerows(rows)

print(f"Wrote {len(rows)} rows to {OUTFILE}")
print(f"\nBottom 10 (worst merging):")
print("\t".join(["Sample", "Pairs", "Joined", "JoinedPct", "NumSourcePairs"]))
for r in rows[:10]:
    print(f"{r['Sample']}\t{r['Pairs']:,}\t{r['Joined']:,}\t{r['JoinedPct']}\t{r['NumSourcePairs']}")
