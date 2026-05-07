#!/bin/bash
#===============================================================================
# File Name    : screen01_mergeR1R2.bash
# Description  : Merge R1 and R2 from Illumina, keep only perfect matches
# Usage        : sbatch screen01_mergeR1R2.bash
# Author       : Aura Ferreiro
# Version      : 2.0
# Created On   : 2024-05-20
# Last Modified: 2026-03-17
#===============================================================================
#SBATCH --array=1-165%20
#SBATCH --time=0-23:00:00
#SBATCH --job-name=perfectmerge
#SBATCH --mem=12G
#SBATCH --cpus-per-task=6
#SBATCH --output=/scratch/rmlab/rmlab_shared3/alferreiro/DMS_PacBio/slurmout/screen/screen01_mergeR1R2_%a.out
#SBATCH --error=/scratch/rmlab/rmlab_shared3/alferreiro/DMS_PacBio/slurmout/screen/screen01_mergeR1R2_%a.err
#SBATCH --mail-type=END
#SBATCH --partition=general

indir="/scratch/rmlab/rmlab_shared3/alferreiro/DMS_PacBio/DantasScreen1_Illumina/CombinedFastq"
outdir="/scratch/rmlab/rmlab_shared3/alferreiro/DMS_PacBio/DantasScreen1_Illumina/R1R2merged"
samplelist="${indir}/sample_basenames.txt"

mkdir -p "${outdir}"

ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${samplelist}")

eval $(/ref/rmlab/software/spack-0.17.2/bin/spack load --sh bbmap@38.63)

echo "Processing ${ID}"

bbmerge.sh in1="${indir}/${ID}_R1.fastq.gz" \
    in2="${indir}/${ID}_R2.fastq.gz" \
    out="${outdir}/${ID}_merged.fastq.gz" \
    outadapter="${outdir}/${ID}_consensus_adaptersequences.txt" \
    verystrict=t \
    threads=${SLURM_CPUS_PER_TASK}

eval $(/ref/rmlab/software/spack-0.17.2/bin/spack unload --sh bbmap@38.63)
