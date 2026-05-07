#!/bin/bash
#===============================================================================
# File Name    : run_dms_pipeline.bash
# Description  : Driver script for the DMS library-QC Snakemake pipeline.
#
# Recommended usage — submit this script itself as a lightweight Slurm "driver"
# job; Snakemake then submits each rule as its own child job:
#
#   cd /scratch/rmlab/rmlab_shared3/alferreiro/DMS_PacBio
#   sbatch shell_scripts/run_dms_pipeline.bash
#
# For a different dataset, pass an alternate config:
#   sbatch shell_scripts/run_dms_pipeline.bash --config config_newrun.yaml
#
# For local execution on an interactive node (small test runs):
#   bash shell_scripts/run_dms_pipeline.bash --mode local
#
# Author       : Aura Ferreiro
# Version      : 1.1
# Created On   : 2026-05-01
# Last Modified: 2026-05-05
#===============================================================================
#
# ── Slurm directives for the driver job (orchestration only, very lightweight) ─
#SBATCH --job-name=snakemake_driver
#SBATCH --time=2-00:00:00          # 2 days — enough for the full pipeline
#SBATCH --mem=2G                   # Snakemake orchestrator uses very little RAM
#SBATCH --cpus-per-task=1
#SBATCH --output=slurmout/snakemake/driver_%j.out
#SBATCH --error=slurmout/snakemake/driver_%j.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ── Parse optional command-line overrides ────────────────────────────────────
CONFIG="config.yaml"
MODE="slurm"
FORCERUN=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)   CONFIG="$2";    shift 2 ;;
        --mode)     MODE="$2";      shift 2 ;;
        --forcerun) FORCERUN="$2";  shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ── Pipeline settings ─────────────────────────────────────────────────────────
SNAKEFILE="Snakefile_library_qc"
CORES=8                   # used only in local mode
SLURM_PARTITION="general"
LOG_DIR="slurmout/snakemake"

# ── Environment setup ─────────────────────────────────────────────────────────
# This block runs on whichever node executes the script (login node when run
# interactively; a compute node when submitted via sbatch).
# It does NOT need to propagate to child jobs — each Snakemake rule sources
# activate_env.sh independently inside its own shell command.
eval "$(/ref/rmlab/software/spack-0.21.1/bin/spack load --sh /g3tt25a)"   # bzip2
eval "$(/ref/rmlab/software/spack-0.21.1/bin/spack load --sh /pscnt54)"   # python 3.11
export LD_LIBRARY_PATH="$SPACK_LD_LIBRARY_PATH"
source "$HOME/alignparse_env/bin/activate"

if ! command -v snakemake &>/dev/null; then
    echo "ERROR: snakemake not found. Install with: pip install snakemake snakemake-executor-plugin-slurm"
    exit 1
fi

mkdir -p "$LOG_DIR"

echo "Driver PID: $$"
echo "Snakemake:  $(snakemake --version)"
echo "Snakefile:  $SNAKEFILE"
echo "Config:     $CONFIG"
echo "Mode:       $MODE"
[[ -n "$FORCERUN" ]] && echo "Force-rerun: $FORCERUN"
echo ""

# ── Dry-run ───────────────────────────────────────────────────────────────────
echo "=== Dry-run ==="
snakemake -s "$SNAKEFILE" --configfile "$CONFIG" --dry-run --quiet 2>&1 | tail -20
echo ""

# ── Optional forcerun flag ────────────────────────────────────────────────────
FORCERUN_FLAG=""
[[ -n "$FORCERUN" ]] && FORCERUN_FLAG="--forcerun $FORCERUN"

# ── Execute ───────────────────────────────────────────────────────────────────
echo "=== Launching pipeline ==="

if [[ "$MODE" == "local" ]]; then
    # Run all steps on the current node (use only on an interactive compute node).
    snakemake -s "$SNAKEFILE" \
        --configfile "$CONFIG" \
        --cores "$CORES" \
        --rerun-incomplete \
        --printshellcmds \
        $FORCERUN_FLAG

elif [[ "$MODE" == "slurm" ]]; then
    # Each Snakemake rule becomes its own sbatch child job.
    snakemake -s "$SNAKEFILE" \
        --configfile "$CONFIG" \
        --executor slurm \
        --jobs 50 \
        --default-resources \
            mem_mb=4000 \
            runtime="4h" \
            slurm_partition="$SLURM_PARTITION" \
        --rerun-incomplete \
        --printshellcmds \
        --latency-wait 60 \
        $FORCERUN_FLAG

else
    echo "ERROR: unknown MODE '$MODE'."
    exit 1
fi

echo ""
echo "Pipeline complete."
