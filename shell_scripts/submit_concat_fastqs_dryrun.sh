#!/bin/bash
#SBATCH --job-name=concat_fastqs_dryrun
#SBATCH --output=/scratch/rmlab/rmlab_shared3/alferreiro/DMS_PacBio/slurmout/screen/concat_fastqs_dryrun_%j.out
#SBATCH --error=/scratch/rmlab/rmlab_shared3/alferreiro/DMS_PacBio/slurmout/screen/concat_fastqs_dryrun_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general

PYTHON=/ref/rmlab/software/alferreiro/micromamba/envs/dms_np126/bin/python
SCRIPT=/scratch/rmlab/rmlab_shared3/alferreiro/DMS_PacBio/python_scripts/build_combined_fastq_manifest.py

echo "Dry-run started: $(date)"
echo "Host: $(hostname)"
echo ""

"$PYTHON" "$SCRIPT" --dry-run

EXIT_CODE=$?
echo ""
echo "Dry-run finished: $(date)  (exit code: $EXIT_CODE)"
exit $EXIT_CODE
