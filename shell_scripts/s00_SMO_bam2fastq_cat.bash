#!/bin/bash
#===============================================================================
# File Name    : s00_SMO_bam2fastq_cat.bash
# Description  : Convert two PacBio BAMs to FASTQ and concatenate into one
#                compressed file. Kinnex (MD969) + Linear (MD990) -> SMO.fastq.gz
# Usage        : sbatch s00_SMO_bam2fastq_cat.bash
# Author       : Aura Ferreiro
# Version      : 1.0
# Created On   : 2026-05-01
#===============================================================================
#SBATCH --time=0-08:00:00
#SBATCH --job-name=SMO_bam2fastq
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8
#SBATCH --output=../slurmout/DBaldridge_SMO/x_SMO_bam2fastq.out
#SBATCH --error=../slurmout/DBaldridge_SMO/y_SMO_bam2fastq.err
#SBATCH --mail-type=END

basedir="/scratch/rmlab/rmlab_shared3/alferreiro/DMS_PacBio/DBaldridge_SMO"
bam1="${basedir}/MD969_Kinnex6Plex/segmented.bam"
bam2="${basedir}/MD990_Linear/m84081_260327_212341_s2.hifi_reads.bc2011.bam"
outfile="${basedir}/SMO.fastq.gz"

eval $(/ref/rmlab/software/spack-0.22/bin/spack load --sh samtools@1.19.2)
eval $(/ref/rmlab/software/spack-0.22/bin/spack load --sh /z3b7vll) # pigz 2.8

# Convert both BAMs to FASTQ and stream into a single compressed output
{ samtools fastq -@ ${SLURM_CPUS_PER_TASK} "${bam1}" ; \
  samtools fastq -@ ${SLURM_CPUS_PER_TASK} "${bam2}" ; } \
  | pigz -p ${SLURM_CPUS_PER_TASK} > "${outfile}"

echo "Read count in SMO.fastq.gz:"
zcat "${outfile}" | echo $((`wc -l`/4))

eval $(/ref/rmlab/software/spack-0.22/bin/spack unload --sh samtools@1.19.2)
eval $(/ref/rmlab/software/spack-0.22/bin/spack unload --sh /z3b7vll)
