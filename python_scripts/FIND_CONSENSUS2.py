#!/usr/bin/env python
#==============================================================================
# File Name :      find_consensus_multi_wofftarget.py
# Description :    Given aligned.csv output of align_parse, find consensus seqs.
# 
# Usage :          python3 find_consensus_multi_wofftarget.py [aligned_output.csv] [outdir] [targetname] [queries]
# Author :         Aura Ferreiro, alferreiro@wustl.edu
# Version :        1.0
# Created On :     2023-12-20
# Last Modified:   2024-02-19
#===============================================================================

# sub diff = none (no limit)
# keep drop

## USAGE DETAILS FOR AP2_alignparse_consensus.bash ##################
# basedir="/path/to/alignparse_out"   #same as outdir in AP1_alignparse_setup.sh
# runname="basename_fq"               #same as run_name in AP1_alignparse_setup.sh, e.g. reads.lenfilt
# libraries="target_seq"              #same as library_name in AP1_alignparse_setup.sh

# indir="${basedir}/${libraries}_align_and_parse/${libraries}_${runname}"
# outdir="${basedir}/${libraries}_consensus_wofftarget"


# mkdir -p ${outdir}

# ID=$(sed -n ${SLURM_ARRAY_TASK_ID}p genes.txt)


# Usage:
# python find_consensus_multi_wofftarget.py <aligned.csv> <outdir> <targetname> <querylist>

# python ../python_scripts/find_consensus_multi_wofftarget.py \
#     --alignment ${indir}/${ID}_PacBio_amplicon_aligned.csv \
#     --outdir ${outdir}/${ID} \
#     --targetname ${ID} \
#     --querylist 'target_seq' 'UMI' 'cPPT' 'mEGFP' \
#     --error_cutoff 1e-3
#####################################################################


# base python
import sys
import os
import warnings
import argparse

# called by plotnine
import numpy
import pandas as pd

# load plotnine and alignparse spack packages in shell script
from plotnine import *

import alignparse.ccs
import alignparse.consensus
import alignparse.minimap2
import alignparse.targets
from alignparse.constants import CBPALETTE

import dms_variants.utils

# Suppress warnings that clutter ouput
warnings.simplefilter("ignore")


# ARGUMENTS ############################################
parser = argparse.ArgumentParser(description='Characterize UMI-barcdoded DMS variants from PacBio CCS data - find consensus sequences.')

parser.add_argument('--alignment', help='Path to output of align_parse.py, a .csv file with minimap2 alignment results')
parser.add_argument('--outdir', help='output directory')
parser.add_argument('--targetname', help='Name of gene target, e.g. MFN2')
parser.add_argument('--querylist', help='Whitespace delimited list of queried elements. Should include "target_seq", "UMI", and any off-target queries for which accuracy and mutations were parsed.',
    nargs="*", type=str, default=["target_seq", "UMI"])
parser.add_argument('--error_cutoff', help='1 - MinAccuracy', type=float, default=1e-3)

# Parse these arguments
args = parser.parse_args()

# Set Outdir
outdir = args.outdir
os.makedirs(outdir, exist_ok=True)



# LOAD DATA ########################################
# load 'aligned' output from align_parse.py 

aligned_path = args.alignment
aligned_df = pd.read_csv(aligned_path, na_filter=None)




# target name 
target = args.targetname

# queries
queries = args.querylist

# error_cutoff 
error_cutoff = args.error_cutoff

# make sure queries for which accuracy is returned are numeric

for query in queries:
    aligned_df[f"{query}_accuracy"] = pd.to_numeric(aligned_df[f"{query}_accuracy"], errors='coerce')


print(f"Finding consensus sequence for {target}.")
print(f"Returning results for the following queries: {queries}")


# detect whether the UMI column is called "UMI" or "barcode"
if "UMI" in queries:
    umi_name = "UMI"
elif "barcode" in queries:
    umi_name = "barcode"
else:
    raise ValueError('`--querylist` must include either "UMI" or "barcode"')





# COMPUTE FLOORED ERROR RATES ######################
# for plotting on log scale
error_rate_floor = 1e-7  # SET error rate floor to this

# Compute errors from reported accuracies by query 
for query in queries:
    accuracy_col = f"{query}_accuracy"
    error_col = f"{query}_error"
    aligned_df = aligned_df.assign(**{
        error_col: lambda x, acc_col=accuracy_col: numpy.clip(1 - x[acc_col], error_rate_floor, None)
        })

# Add library column
aligned_df['library'] = target


# EXCLUDE BARCODES/GENES WITH ERROR RATES > CUTOFF##
#error_cutoff = 1e-3


error_columns = [f"{query}_error" for query in queries]
aligned_errordist_df = aligned_df.melt(
    id_vars=["library"],
    value_vars=error_columns,
    var_name="feature_type",
    value_name="error rate"
)

errordist_path = f"{outdir}/errordistribution_alignedCCSs.csv"
aligned_errordist_df.to_csv(errordist_path, index=False)

library_count = aligned_errordist_df['library'].nunique()
feature_count = aligned_errordist_df['feature_type'].nunique()

# Plot distribution of error rates
p_err = (
    ggplot(
        aligned_errordist_df,
        aes("error rate"),
    )
    + geom_histogram(bins=25)
    + geom_vline(xintercept=error_cutoff, linetype="dashed", color=CBPALETTE[1])
    + facet_grid("library ~ feature_type")
    + theme(figure_size=(3 * feature_count, 5 * library_count),
        axis_text=element_text(size=12),
        axis_title=element_text(size=15))
    + ylab("number of CCSs")
    + scale_x_log10()
)



plotpath = f"{outdir}/error_rate_distribution.png"
ggsave(p_err, filename=plotpath, limitsize=False)



# Mark the CCSs that pass the error filter (based just on target seq and UMI queries)
aligned_df = aligned_df.assign(
    retained=lambda x: (
        (x["target_seq_error"] <= error_cutoff) &
        (x[f"{umi_name}_error"]   <= error_cutoff)
    )
)



# Summarize the CCSs that pass the error filter
alignpass_df = (
    aligned_df.groupby(["library", "retained"])
    .size()
    .rename("number of CCSs")
    .reset_index()
)

alignpass_path = f"{outdir}/Num_CCSs_passing_errorfilter.csv"
alignpass_df.to_csv(alignpass_path, index=False)




# Plot number of CCSs per barcode (for retained CCSs)
# MFN2 has long tail
max_CCSs = 100  # in plot, group barcodes with >= this many CCSs
library_count = aligned_df['library'].nunique()

p_CCSsperBC = (
    ggplot(
        aligned_df.query("retained")
        .groupby(["library", f"{umi_name}_sequence"])
        .size()
        .rename("nseqs")
        .reset_index()
        .assign(nseqs=lambda x: numpy.clip(x["nseqs"], None, max_CCSs)),
        aes("nseqs"),
    )
    + geom_bar()
    + facet_wrap("~ library", nrow=1)
    + theme(
        figure_size=(5 * library_count, 5),
        panel_grid_major_x=element_blank(),
        axis_text=element_text(size=12),
        axis_title=element_text(size=15)
    )
    + ylab("number of barcodes")
    + xlab("number of CCSs per barcode")
)

CCSsperBC_path = f"{outdir}/CCSs_per_barcode.png"
ggsave(p_CCSsperBC, filename=CCSsperBC_path)



# ESTIMATE EMPIRICAL ACCURACY BY COMPARING CCSs WITH SAME UMI #######


# Add information about whether the target seqs have indels
aligned_df = alignparse.consensus.add_mut_info_cols(
    aligned_df, mutation_col="target_seq_mutations", n_indel_col="n_indels"
)


aligned_df = aligned_df.assign(has_indel=lambda x: x["n_indels"] > 0)


# Summarize +/- indel numbers among retained CCSs and save
numCCSs_wIndels = (
    aligned_df.query("retained")
    .groupby(["library", "has_indel"])
    .size()
    .rename("number_CCSs")
    .reset_index()
)


numindels_path = f"{outdir}/NumCCSs_wIndels.csv"
numCCSs_wIndels.to_csv(numindels_path, index=False)

# WITH INDELS
# and print to standard out

#print("Empirical Accuracy With Indels")

#w_indels = alignparse.consensus.empirical_accuracy(
#    aligned_df.query("retained"), mutation_col="target_seq_mutations",
#    accuracy_col='target_seq_empirical_accuracy', group_cols=f"{umi_name}_sequence"
#)

#print(w_indels)

# WITHOUT INDELS
# and print to standard out

#print("Empirical Accuracy Without Indels")

#wo_indels = alignparse.consensus.empirical_accuracy(
#    aligned_df.query("retained & not has_indel"), mutation_col="target_seq_mutations",
#    accuracy_col='target_seq_empirical_accuracy', group_cols=f"{umi_name}_sequence"
#)

#print(wo_indels)



# BUILD THE CONSENSUS SEQUENCES ####################
mutqueries = [query for query in queries if query != umi_name]

for query in mutqueries:

    #consensus, dropped = alignparse.consensus.simple_mutconsensus(
    #    aligned_df.query("retained"),
    #    group_cols=("library", "UMI_sequence"),
    #    mutation_col=f"{query}_mutations",
    #)

    consensus, dropped = alignparse.consensus.simple_mutconsensus(
        aligned_df.query("retained"),
        group_cols=("library", f"{umi_name}_sequence"),
        mutation_col=f"{query}_mutations",
        max_sub_diffs=None,
    )

    # Add info on which mutations are just substitutions and which have indels
    consensus = alignparse.consensus.add_mut_info_cols(
        consensus,
        mutation_col=f"{query}_mutations",
        sub_str_col="substitutions",
        n_indel_col="number_of_indels",
        overwrite_cols=True,
    )

    # save 
    consensus_path = f"{outdir}/CONSENSUS_mutations_{query}.csv"
    consensus.to_csv(consensus_path, index=False)
 
    # save dropped too
    dropped_path = f"{outdir}/CONSENSUS_mutations_{query}_dropped.csv"
    dropped.to_csv(dropped_path, index=False)


    # Summarize the barcodes that were dropped and save
    dropped_summ = (dropped.groupby("drop_reason").size().rename("number_barcodes").reset_index())

    droppedsumm_path = f"{outdir}/CONSENSUS_droppedbarcodes_{query}.csv"
    dropped_summ.to_csv(droppedsumm_path, index=False)



# CONSENSUS, FILTERED TO REMOVE INDELS, CAN BE FURTHER PROCESSED
# BY NEXT SCRIPT: VARIANT_ANALYSIS.PY





