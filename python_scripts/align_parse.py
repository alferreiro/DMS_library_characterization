#!/usr/bin/env python
#==============================================================================
# File Name :      align_parse.py
# Description :    Extract UMIs, count/summarize, using alignparse package.
#                  Part I of the alignparse pipeline.
# Usage :          python3 align_parse.py <arguments>
# Author :         Aura Ferreiro, alferreiro@wustl.edu
# Version :        1.2
# Created On :     2023-12-11
# Last Modified:   2026-05-07
#   v1.1: added --fastq full-path arg, .fastq.gz support
#   v1.2: --minimap2 CLI arg (was hardcoded); targeted warning suppression;
#         required=True for all mandatory args
#===============================================================================

import sys
import os
import warnings
import argparse

import numpy
import pandas as pd
from plotnine import *

import alignparse.ccs
import alignparse.consensus
import alignparse.minimap2
import alignparse.targets
from alignparse.constants import CBPALETTE

import dms_variants.utils

# Suppress DeprecationWarning and FutureWarning from third-party libraries
# (alignparse, dms_variants, plotnine emit many of these internally).
# RuntimeWarning and UserWarning remain visible as they may indicate real problems.
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# ARGUMENTS ####################################################################
parser = argparse.ArgumentParser(
    description='Characterize UMI-barcoded DMS variants from PacBio CCS data.')

parser.add_argument('--outdir',       required=True,  help='output directory')
parser.add_argument('--target_gb',    required=True,  help='annotated genbank file')
parser.add_argument('--library_name', required=True,  help='library label, e.g. "target_seq"')
parser.add_argument('--yaml',         required=True,  help='alignparse config yaml')
parser.add_argument('--minimap2',     required=True,
                    help='path to minimap2 executable')

# Two ways to specify the fastq: full path OR run_name + fq_dir
fastq_group = parser.add_mutually_exclusive_group(required=True)
fastq_group.add_argument('--fastq',    help='full path to fastq.gz file')
fastq_group.add_argument('--run_name', help='fastq basename (without extension); requires --fq_dir')
parser.add_argument('--fq_dir', help='directory containing the fastq file (used with --run_name)')

args = parser.parse_args()

# Resolve fastq path and run_name
if args.fastq:
    fastq_path = args.fastq
    run_name = os.path.basename(fastq_path)
    for ext in ['.fastq.gz', '.fq.gz', '.fastq', '.fq']:
        if run_name.endswith(ext):
            run_name = run_name[:-len(ext)]
            break
else:
    if not args.fq_dir:
        parser.error("--fq_dir is required when using --run_name")
    run_name = args.run_name
    for ext in ['.fq.gz', '.fastq.gz']:
        fastq_path = os.path.join(args.fq_dir, run_name + ext)
        if os.path.exists(fastq_path):
            break
    else:
        raise FileNotFoundError(
            f"Could not find {run_name}.fq.gz or {run_name}.fastq.gz in {args.fq_dir}")

outdir         = args.outdir
target_genbank = args.target_gb
libraries      = args.library_name
config_file    = args.yaml

os.makedirs(outdir, exist_ok=True)

pacbio_runs = pd.DataFrame({
    "name":    [run_name],
    "library": [libraries],
    "fastq":   [fastq_path],
})

# QC summary plots
ccs_summaries = alignparse.ccs.Summaries(pacbio_runs, report_col=None)
for stat in ["length", "passes", "accuracy"]:
    if ccs_summaries.has_stat(stat):
        p = ccs_summaries.plot_ccs_stats(stat)
        p = p + theme(panel_grid_major_x=element_blank(),
                      axis_text=element_text(size=12),
                      axis_title=element_text(size=15))
        ggsave(p, filename=f"{outdir}/{stat}.summary.png", width=10, height=10)
        ccs_summaries.ccs_stats(stat).to_csv(f"{outdir}/{stat}.csv", index=False)
    else:
        print(f"No {stat} statistics available.")

# Build targets object
targets = alignparse.targets.Targets(
    seqsfile=target_genbank, feature_parse_specs=config_file)

genbank_plot = targets.plot(ax_width=10)
genbank_plot.savefig(f"{outdir}/genbank_targetmap.png")

print(tuple(pacbio_runs["name"]))
print('Parsing configuration:\n')
print(targets.feature_parse_specs("yaml"))

# minimap2 mapper — path provided via --minimap2 argument
mapper = alignparse.minimap2.Mapper(
    alignparse.minimap2.OPTIONS_CODON_DMS,
    prog=args.minimap2,
    retain_tags=None)
print(f"minimap2 {mapper.version}: " + " ".join(mapper.options))

# Align and parse
align_and_parse_outdir = os.path.join(outdir, f"{libraries}_align_and_parse")
readstats, aligned, filtered = targets.align_and_parse(
    df=pacbio_runs,
    mapper=mapper,
    outdir=align_and_parse_outdir,
    name_col="name",
    queryfile_col="fastq",
    group_cols=["library"],
    to_csv=True,
    overwrite=True,
    ncpus=-1,
)
