#!/bin/bash
# activate_env.sh — Source this in every SLURM job

# Load the same Spack environment as in setup
eval "$(/ref/rmlab/software/spack-0.21.1/bin/spack load --sh /g3tt25a)"   # bzip2 built with gcc@8.5.0
eval "$(/ref/rmlab/software/spack-0.21.1/bin/spack load --sh /pscnt54)"   #python 3.11 built with gcc@8.5.0

# Set required paths
export LD_LIBRARY_PATH="$SPACK_LD_LIBRARY_PATH"
export C_INCLUDE_PATH="$SPACK_C_INCLUDE_PATH"

# Activate virtualenv
source "$HOME/alignparse_env/bin/activate"