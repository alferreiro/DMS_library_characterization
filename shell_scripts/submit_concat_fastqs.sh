#!/bin/bash
#SBATCH --job-name=concat_fastqs
#SBATCH --output=/scratch/rmlab/rmlab_shared3/alferreiro/DMS_PacBio/slurmout/screen/concat_fastqs_%j.out
#SBATCH --error=/scratch/rmlab/rmlab_shared3/alferreiro/DMS_PacBio/slurmout/screen/concat_fastqs_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --partition=general

PYTHON=/ref/rmlab/software/alferreiro/micromamba/envs/dms_np126/bin/python
SCRIPT=/scratch/rmlab/rmlab_shared3/alferreiro/DMS_PacBio/python_scripts/build_combined_fastq_manifest.py

echo "Job started: $(date)"
echo "Host: $(hostname)"
echo ""

"$PYTHON" "$SCRIPT"

EXIT_CODE=$?
echo ""
echo "Job finished: $(date)  (exit code: $EXIT_CODE)"
exit $EXIT_CODE
