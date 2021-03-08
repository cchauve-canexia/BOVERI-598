#!/usr/bin/env python3

"""
Check the presence of a set of indels in observed indels
"""
import argparse
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import yaml

# Dataframes column keys: should be imported from vcf_utils.py
from analyze_variants_utils import (
    SAMPLE,
    RUN_ID,
    VAF,
    EXP_VAF,
    CHR,
    POS,
    REF,
    ALT,
    NG,
    W_SCORE,
    SCORE,
    COMPLEXITY,
    SUPPORT,
    OVERLAP,
    CONTROL,
)

# Used columns of blacklist files
INDEL_ID = 'id'

# ------------------------------------------------------------------------------
# Auxiliary functions


def coord_to_del(chr, start, end, manifest_df):
    """
    Given an interval chr:start-end, returns the corresponding sequence
    extended by one base to the left
    """
    for _, amplicon in manifest_df.iterrows():
        test_1 = amplicon['Chr'] == chr
        test_2 = int(amplicon['Start']) < start
        test_3 = int(amplicon['End']) >= end
        if test_1 and test_2 and test_3:
            pos = start - int(amplicon['Start']) - 1
            ref_len = end - start + 1
            ref = amplicon['Amplicon'][pos:pos + ref_len + 2]
            alt = amplicon['Amplicon'][pos:pos + 1]
            return {CHR: chr, POS: pos, REF: ref, ALT: alt}

def read_indels(indels_tsv_file, manifest_df, sep='\t'):
    """
    Extracts from a CSV file a list of (chr, pos, ref, alt) describing
    indels
    :param: blacklist_csv_file (str): path to a CSV file describing indels in
    the format chr,pos,ref,alt
    :param: manifest_df (DataFrame): amplicon manifest
    :param: sep (str): separator in CSV file
    :return: list(dict(str, str, int, str, str)): indexed by
    CHR, POS, REF, ALT
    """
    indels = []
    indels_df = pd.read_csv(indels_tsv_file, sep='\t')
    for _, row in indels_df.iterrows():
        indel_str = row[INDEL].split(':')
        indel = {
            CHR: indel_str[0]
        }
        if len(indel_str) == 4:
            indel[POS] = int(indel_str[1])
            indel[REF], indel[ALT] = indel_str[2], indel_str[3]
        else:
            assert indel_str[2] == 'del'
            [start, end] = [int(x) for x in indel_str[1].split(',')]
            indel_dict = coord_to_del(indel_str[0], start, end, manifest_df)
            for key in [POS, REF, ALT]:
                indel[key] = indel_dict[key]
        indels.append(indel)
    return indels

def get_runs_data(run_id_list):
    """
    Returns the dataframes of observed indels for run in run_id_list and
    sample from sample_list from
    results files
    :param: run_id_list (list(str)): ID of the runs to consider
    :return DataFrame: dataframe of observed indels
    """
    observed_indels_df_list = []
    for run_id in run_id_list:
        # Reading the pipeline results
        run_indels_file = os.path.join(
            'results', run_id, f"{run_id}_indels.tsv"
        )
        observed_indels_df = pd.read_csv(run_indels_file, sep='\t')
        # Reformatting the sample column to match the true indels format
        observed_indels_df[SAMPLE] = observed_indels_df.apply(
            lambda row: row[SAMPLE].split('_S')[0], axis=1
        )
        observed_indels_df_1 = observed_indels_df.loc[
            ~observed_indels_df[SAMPLE].str.startswith('B')
        ].round(3)
        observed_indels_df_list.append(observed_indels_df_1)
        # Excluding indels in the blasklist
    all_observed_indels_df = pd.concat(observed_indels_df_list)
    return all_observed_indels_df


# ------------------------------------------------------------------------------
if __name__ == "__main__":
    # Analysis parameters
    # Input file
    ARGS_RUNS_FILE = ['runs_file', None, 'Runs CSV file']
    ARGS_TSV_INDELS_FILE = ['tsv_indels_file', None, 'Checked indels file']
    parser = argparse.ArgumentParser(description='Indels testing: checking indels occurrences')
    parser.add_argument(ARGS_RUNS_FILE[0],
                        type=str,
                        help=ARGS_RUNS_FILE[2])
    parser.add_argument(ARGS_TSV_INDELS_FILE[0],
                        type=str,
                        help=ARGS_TSV_INDELS_FILE[2])
    args = parser.parse_args()
    out_pref = args.runs_file.replace('data', 'results').replace('.tsv', '')
    # Reading runs
    RUNS_DF = pd.read_csv(args.runs_file, sep=',', header=None, names=['run', RUN_ID])
    # Reading checked indels file
    CHECKED_INDELS_DF = pd.read_csv(args.tsv_indels_file, sep='\t')
    # List of runs and samples with at least one expected indel
    RUNS_LIST = list(RUNS_DF[RUN_ID].unique())
    RUNS_INDELS_DF = get_runs_data(RUNS_LIST)
    SAMPLES_LIST = list(
        RUNS_INDELS_DF.loc[
            ~RUNS_INDELS_DF[SAMPLE].str.startswith('B')
        ][SAMPLE].unique()
    )
    NB_SAMPLES = len(SAMPLES_LIST)

    # Processing checked indels
    out_file = open(f"{out_pref}_out_1.tsv", 'w')
    min_nb_occurrences = NB_SAMPLES
    min_mean_vaf = 1.0
    checked_indels_list = []
    out_file.write(f"# Number of non-Blank samples: {len(SAMPLES_LIST)}")
    out_file.write('\nindel\tnb_occ\tmin_vaf\tmax_vaf\tmean_vaf\tmin_ctrl\tmax_ctrl\tmean_ctrl')
    for _, indel_row in CHECKED_INDELS_DF.iterrows():
        indel_1 = indel_row[INDEL_ID].split(':')
        indel = {
            CHR: indel_1[0],
            POS: int(indel_1[1]),
            REF: indel_1[2],
            ALT: indel_1[3]
        }
        checked_indels_list.append((indel[CHR], indel[POS], indel[REF], indel[ALT]))
        indel_df = RUNS_INDELS_DF.loc[
            (RUNS_INDELS_DF[CHR]==indel[CHR]) &
            (RUNS_INDELS_DF[POS]==indel[POS]) &
            (RUNS_INDELS_DF[REF]==indel[REF]) &
            (RUNS_INDELS_DF[ALT]==indel[ALT])
        ]
        nb_occurrences = len(indel_df.index)
        if nb_occurrences < min_nb_occurrences:
            min_nb_occurrences = nb_occurrences
        max_vaf = round(np.max(indel_df[VAF]), 2)
        min_vaf = round(np.nanmin(indel_df[VAF]), 2)
        mean_vaf = round(np.mean(indel_df[VAF]), 2)
        if mean_vaf < min_mean_vaf:
            min_mean_vaf = mean_vaf
        max_control = round(np.max(indel_df[CONTROL]), 2)
        min_control = round(np.nanmin(indel_df[CONTROL]), 2)
        mean_control = round(np.mean(indel_df[CONTROL]), 2)
        out_file.write(f"\n{indel[CHR]}:{indel[POS]}:{indel[REF]}:{indel[ALT]}\t{nb_occurrences}\t{min_vaf}\t{max_vaf}\t{mean_vaf}\t{min_control}\t{max_control}\t{mean_control}")
    out_file.close()

    # Statistics of indels with similar features
    out_file = open(f"{out_pref}_out_2.tsv", 'w')
    out_file.write('indel\tnb_occ\tmin_vaf\tmax_vaf\tmean_vaf\tmin_ctrl\tmax_ctrl\tmean_ctrl')
    INDELS_DF = RUNS_INDELS_DF.copy()
    INDELS_GROUPS = INDELS_DF.groupby(by=[CHR, POS, REF, ALT])
    for indel, indel_df in INDELS_GROUPS:
        if indel not in checked_indels_list:
            nb_indels = len(list(indel_df.index))
            mean_vaf = np.mean(indel_df[VAF])
            if nb_indels >= min_nb_occurrences and mean_vaf >= min_mean_vaf:
                max_vaf = round(np.max(indel_df[VAF]), 2)
                min_vaf = round(np.nanmin(indel_df[VAF]), 2)
                mean_vaf = round(np.mean(indel_df[VAF]), 2)
                max_control = round(np.max(indel_df[CONTROL]), 2)
                min_control = round(np.nanmin(indel_df[CONTROL]), 2)
                mean_control = round(np.mean(indel_df[CONTROL]), 2)
                out_file.write(f"{indel[0]}:{indel[1]}:{indel[2]}:{indel[3]}\t{nb_occurrences}\t{min_vaf}\t{max_vaf}\t{mean_vaf}\t{min_control}\t{max_control}\t{mean_control}")
    out_file.close()

    # Generating a scatter plot for all indels occurrences versus mean VAF
    # and mean control penalty, recording multiplicity of all indels
    NB_OCC = defaultdict(int)
    NB_COL, VAF_COL, CTRL_COL, STATUS_COL = 'nb', 'vaf', 'ctrl', 'status'
    SUPP_COL, OV_COL, COMP_COL = 'support', 'overlap', 'complexity'
    INDELS_GROUPS_DICT = {}
    for indel, indel_df in INDELS_GROUPS:
        nb_indels = len(list(indel_df.index))
        mean_vaf = indel_df[VAF].mean()
        mean_ctrl = indel_df[CONTROL].mean()
        mean_supp = indel_df[SUPPORT].mean()
        mean_ov = indel_df[OVERLAP].mean()
        mean_comp = indel_df[COMPLEXITY].mean()
        if indel in checked_indels_list:
            status = 'blacklist'
        else:
            status = 'non-blakclist'
        INDELS_GROUPS_DICT[indel] = {
            NB_COL: nb_indels, VAF_COL: mean_vaf,
            CTRL_COL: mean_ctrl, STATUS_COL: status,
            SUPP_COL: mean_supp, OV_COL: mean_ov, COMP_COL: mean_comp
        }
        NB_OCC[nb_indels] += 1
    OCC_KEYS = list(NB_OCC.keys())
    OCC_KEYS.sort()
    out_file = open(f"{out_pref}_out_3.tsv", 'w')
    out_file.write('nb_occurrences\tnb_indels')
    for nb_occ in OCC_KEYS:
        out_file.write(f"\n{nb_occ}\t{NB_OCC[nb_occ]}")
    out_file.close()

    INDELS_GROUPS_DF = pd.DataFrame.from_dict(
        INDELS_GROUPS_DICT,
        orient='index',
        columns=[NB_COL, VAF_COL, CTRL_COL, STATUS_COL, SUPP_COL, OV_COL, COMP_COL]
    )
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(20, 20))
    vaf_plot = sns.scatterplot(
        data=INDELS_GROUPS_DF,
        x=NB_COL, y=VAF_COL, hue=STATUS_COL,
        ax=axes[0], alpha=0.5
    )
    vaf_plot = sns.scatterplot(
        data=INDELS_GROUPS_DF,
        x=NB_COL, y=CTRL_COL, hue=STATUS_COL,
        ax=axes[1], alpha=0.5
    )
    axes[0].set_title(f"All indels, mean VAF versus nb occ.")
    axes[1].set_title(f"All indels, mean control penalty versus nb occ.")
    plt.savefig(f"{out_pref}_out_4.png")

    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(20, 30))
    vaf_plot = sns.scatterplot(
        data=INDELS_GROUPS_DF,
        x=NB_COL, y=SUPP_COL, hue=STATUS_COL,
        ax=axes[0], alpha=0.5
    )
    vaf_plot = sns.scatterplot(
        data=INDELS_GROUPS_DF,
        x=NB_COL, y=OV_COL, hue=STATUS_COL,
        ax=axes[1], alpha=0.5
    )
    vaf_plot = sns.scatterplot(
        data=INDELS_GROUPS_DF,
        x=NB_COL, y=COMP_COL, hue=STATUS_COL,
        ax=axes[2], alpha=0.5
    )
    axes[0].set_title(f"All indels, mean support penalty versus nb occ.")
    axes[1].set_title(f"All indels, mean overlap penalty versus nb occ.")
    axes[2].set_title(f"All indels, mean complexity penalty versus nb occ.")
    plt.savefig(f"{out_pref}_out_5.png")
