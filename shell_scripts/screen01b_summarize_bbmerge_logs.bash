#!/bin/bash
#===============================================================================
# File Name    : screen01b_summarize_bbmerge_logs.bash
# Description  : Count read pairs and joined read pairs after bbmerge, from logs,
#                summing by sample
# Usage        : sbatch screen01b_summarize_bbmerge_logs.bash
# Author       : Aura Ferreiro
# Version      : 2.0
# Created On   : 2025-09-10
# Last Modified: 2026-03-17
#===============================================================================
#SBATCH --array=1-1
#SBATCH --time=0-01:00:00
#SBATCH --job-name=summarize_bbmerge
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/rmlab/rmlab_shared3/alferreiro/DMS_PacBio/slurmout/screen/screen01b_summarizebbmerge_%a.out
#SBATCH --error=/scratch/rmlab/rmlab_shared3/alferreiro/DMS_PacBio/slurmout/screen/screen01b_summarizebbmerge_%a.err
#SBATCH --mail-type=END
#SBATCH --partition=general

set -euo pipefail

LOGDIR="/scratch/rmlab/rmlab_shared3/alferreiro/DMS_PacBio/slurmout/screen"
PATTERN="screen01_mergeR1R2_*.err"
OUTDIR="/scratch/rmlab/rmlab_shared3/alferreiro/DMS_PacBio/DantasScreen1_Illumina/R1R2merged"

PERFILE="${OUTDIR}/bbmerge_perfile.tsv"
GROUPED="${OUTDIR}/bbmerge_grouped.tsv"

printf "LogFile\tSample\tPairs\tJoined\tJoinedPct\n" > "$PERFILE"

shopt -s nullglob
for f in "$LOGDIR"/$PATTERN; do
  # Extract sample name from the in1= argument on the first line of the log
  sample="$(grep -oE 'in1=[^ ]+' "$f" | head -n1 | sed 's|.*/||; s/_R1\.fastq\.gz$//')"

  if [[ -z "${sample:-}" ]]; then
    printf "WARNING: could not parse sample name in %s\n" "$f" >&2
    continue
  fi

  # Numbers appear as "Pairs:\t#######" or with commas/percentages; strip accordingly
  pairs_line="$(grep -m1 -E '^[[:space:]]*Pairs' "$f" || true)"
  joined_line="$(grep -m1 -E '^[[:space:]]*Joined' "$f" || true)"
  pairs="$(printf "%s" "$pairs_line"  | grep -oE '[0-9][0-9,]*' | head -n1 | tr -d ',' || true)"
  joined="$(printf "%s" "$joined_line" | grep -oE '[0-9][0-9,]*' | head -n1 | tr -d ',' || true)"
  joined_pct="$(printf "%s" "$joined_line" | grep -oE '[0-9]+\.[0-9]+%' | head -n1 || true)"

  if [[ -z "${pairs:-}" || -z "${joined:-}" ]]; then
    printf "WARNING: could not parse Pairs/Joined fields in %s\n" "$f" >&2
    continue
  fi

  printf "%s\t%s\t%s\t%s\t%s\n" \
    "$(basename "$f")" "$sample" "$pairs" "$joined" "${joined_pct:-NA}" >> "$PERFILE"
done

n_parsed=$(( $(wc -l < "$PERFILE") - 1 ))
printf "Parsed %d log files.\n" "$n_parsed"

# Sort by sample name and summarise (each sample has one log file, but grouping
# handles the edge case where reruns produced duplicate entries)
awk -F'\t' 'NR>1 {pairs[$2]+=$3; joined[$2]+=$4} END {
  print "Sample\tPairs\tJoined\tJoinedPct";
  n=asorti(pairs, idx);
  for (i=1; i<=n; i++) {
    s=idx[i];
    pct=(pairs[s]>0) ? sprintf("%.3f%%", joined[s]/pairs[s]*100) : "NA";
    printf "%s\t%d\t%d\t%s\n", s, pairs[s], joined[s], pct;
  }
}' "$PERFILE" > "$GROUPED"

echo "Wrote per-file summary : $PERFILE"
echo "Wrote grouped summary  : $GROUPED"
