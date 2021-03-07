#!/usr/bin/env python3

"""
Check the presence of a set of indels in observed indels
"""
import argparse
import numpy as np
import os
import pandas as pd
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
    ARGS_MIN_VAF = ['-m', '--min_vaf', 'Min VAF to consider indels']
    parser = argparse.ArgumentParser(description='Indels testing: checking indels occurrences')
    parser.add_argument(ARGS_RUNS_FILE[0],
                        type=str,
                        help=ARGS_RUNS_FILE[2])
    parser.add_argument(ARGS_TSV_INDELS_FILE[0],
                        type=str,
                        help=ARGS_TSV_INDELS_FILE[2])
    parser.add_argument(ARGS_MIN_VAF[0],
                        ARGS_MIN_VAF[1],
                        type=float,
                        help=ARGS_MIN_VAF[2])
    args = parser.parse_args()
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
    min_nb_occurrences = NB_SAMPLES
    min_mean_vaf = 1.0
    # Processing indels for all parameters and threshold settings
    checked_indels_list = []
    print(f"Number of non-Blank samples: {len(SAMPLES_LIST)}")
    print('indel\tnb_occ\tmin_vaf\tmax_vaf\tmean_vaf\tmin_ctrl\tmax_ctrl\tmean_ctrl')
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
        print(f"{indel[CHR]}:{indel[POS]}:{indel[REF]}:{indel[ALT]}\t{nb_occurrences}\t{min_vaf}\t{max_vaf}\t{mean_vaf}\t{min_control}\t{max_control}\t{mean_control}")

    # Looking for potential widespread indels
    if args.min_vaf is not None:
        print('\nindel\tnb_occ\tmin_vaf\tmax_vaf\tmean_vaf\tmin_ctrl\tmax_ctrl\tmean_ctrl')
        INDELS_DF = RUNS_INDELS_DF.loc[RUNS_INDELS_DF[VAF]>=args.min_vaf]
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
                    print(f"{indel[0]}:{indel[1]}:{indel[2]}:{indel[3]}\t{nb_occurrences}\t{min_vaf}\t{max_vaf}\t{mean_vaf}\t{min_control}\t{max_control}\t{mean_control}")
