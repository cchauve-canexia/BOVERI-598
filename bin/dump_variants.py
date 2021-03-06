#!/usr/bin/env python3

"""
Filtering part of the pipeline
"""
import argparse
import glob
import os
import pandas as pd
import time

from bin.alignments_utils import AlignmentsSet
from bin.amplicons_utils import AmpliconsManifest
from bin.clusters_utils import ClusterReadsSet
from bin.features_utils import (
    SCORE,
)
from bin.parameters_utils import (
    CONTROL_BLANK,
    CONTROL_NF,
    Parameters,
)
from bin.pipeline_main import (
    INDELS_VCF_EXT,
)
from bin.pipeline_file_names import (
    clusters_txt_file,
    alignments_txt_file,
    variants_graph_txt_file,
    vcf_sample_file,
)
from bin.postprocessing_utils import (
    aggregate_variants,
    filter_min_quality,
    filter_nb_gaps,
)
from postprocessing_utils import (
    add_confidence_score,
    add_control_score,
    add_support_score,
    add_complexity_score,
    add_overlap_score,
    get_control_samples_feature,
    SCORE_DEFAULT,
)
from bin.variants_graph_utils import VariantsGraph
from bin.variants_utils import (
    INDEL,
    SNP,
)
from vcf_utils import (
    dump_vcf_to_tsv,
    vcf_header,
    run_snpeff,
    vcf_write_variants_features,
    vcf_write_df,
    vcf_import_df,
    CONTROL,
    COMPLEXITY,
    OVERLAP,
    SUPPORT,
    VCF_SORT_POS,
    SAMPLE,
    CHR_COL,
    POS_COL,
    INFO_COL,
    ID_COL,
    RUN_ID,
    RUN_NAME,
)


# Parameters for penalty score functions
COMP_KMIN = 2 # Min k-mers considered for sequence complexity score
COMP_KMAX = 3 # Max k-mers considered for sequence complexity score
COMP_L = 5 # Number of flanking bases for sequence complexity score

# Extension for VCF files with scores computed
FILTERED_INDELS_VCF_EXT = '_indels_filtered'
MANIFEST_V4_FILE = 'CG001v4.0_Amplicon_Manifest_Panel4.0.3_20181101.tsv'

def apply_amplicon_filters(parameters, sample_id, amplicon_id, manifest):
    """
    Filters on all calls from a single amplicon on the variants graph
    :param: parameters (Parameters): pipeline parameters
    :param: in_list (list(Variant, VariantFeatures)): input list of variants
    and their features
    :param: sample_id (str): sample ID
    :param: amplicon ID (str): amplicon ID
    :param: manifest (AmpliconsManifest): amplicon manifest
    :return:  list(Variant, VariantFeatures): list of un-filtered indels
    and their features
    """
    # Loading data needed for filters, amplicon sequence, clusters,
    # alignments, variants-graph
    amplicon_seq = manifest.get_amplicon(amplicon_id).get_seq()
    clusters_file = clusters_txt_file(parameters, sample_id, amplicon_id)
    clusters = ClusterReadsSet.from_file(None, clusters_file)
    alignments_file = alignments_txt_file(
        parameters, sample_id, amplicon_id
    )
    alignments = AlignmentsSet.from_file(None, alignments_file)
    graph_file = variants_graph_txt_file(parameters, sample_id, amplicon_id)
    graph = VariantsGraph.from_file(None, graph_file)
    graph.set_variants_features(amplicon_seq, amplicon_id, clusters)
    # Applying the filters to the graph:
    # eliminating edges from alignments  with too many gaps
    filter_nb_gaps(parameters, clusters, alignments, graph)
    filter_min_quality(parameters, clusters, alignments, graph, sample_id)
    # Reset features to update coverage
    graph.set_variants_features(amplicon_seq, amplicon_id, clusters)
    graph.set_variants_vaf()
    graph.set_chr_coordinates(manifest)
    new_variants_list = graph.get_variants_and_features()
    # Creating the lists of un-filtered indels and snps
    amp_indels_list = []
    for (variant, features) in new_variants_list:
        if  variant.get_type() in INDEL:
            amp_indels_list.append((variant, features))
    return amp_indels_list

def process_sample(
    parameters, sample_id, run_id, run_name, amplicon_list,
    manifest, control_calls, out_vcf_file, control_sample=True
):
    """
    Filter SNP and indels calls for a sample, using control_calls to filter
    with control_samples if control_sample is True and populating
    control_calls otherwise.
    Export the resulting variants in two VCF files, one prior to run snpEff
    and one after running snpEff.

    :param: parameters (Parameters): pipeline parameters
    :param: sample_id (str): ID of the sample
    :param: amplicon_list (list(str)): list of considered amplicons
    :param: manifest (AmpliconsManifest): amplicons manifest
    :param: control_calls (dict(sample_id ->
    list((Variant, VariantFeatures)))): control variants
    :param: control_sample (bool): True if sample_id is a control sample
    """
    def reformat_info(row):
        """
        Reformat the INFO field to integrate the actual scores recorded
        in other fields
        """
        info_list = row[INFO_COL].split(';')
        info_dict = {x.split('=')[0]: x.split('=')[1] for x in info_list}
        info_dict[SCORE] = row[SCORE]
        info_dict[COMPLEXITY] = row[COMPLEXITY]
        info_dict[SUPPORT] = row[SUPPORT]
        info_dict[OVERLAP] = row[OVERLAP]
        info_dict[CONTROL] = row[CONTROL]
        info_str = ';'.join([f"{k}={v}" for k, v in info_dict.items()])
        return info_str

    # Creating lists of all indels from the variants graph
    # filtered to remove support by alignments with a number of gaps
    # above filter_parameters[FILTER_NB_GAPS]
    sample = sample_id.split('_S')[0]
    in_indels_list = []
    for amplicon_id in amplicon_list:
        amp_indels_list = apply_amplicon_filters(
            parameters, sample_id, amplicon_id, manifest
        )
        in_indels_list += amp_indels_list
    # Filtering indels per sample
    aggregated_indels_list_1 = aggregate_variants(in_indels_list)
    # Extra features default values
    score_dict = {}
    for key in [SCORE, COMPLEXITY, SUPPORT, OVERLAP, CONTROL]:
        score_dict[key] = SCORE_DEFAULT
    score_dict[SAMPLE] = sample
    score_dict[RUN_ID] = run_id
    score_dict[RUN_NAME] = run_name
    if control_sample:
        # Recording but not writing indel calls
        aggregated_indels_list = [
            (variant, features, score_dict.copy())
            for (variant, features) in aggregated_indels_list_1
        ]
        control_calls[sample_id] = aggregated_indels_list
        out_indels_ext = f"{FILTERED_INDELS_VCF_EXT}_aux"
        out_indels_vcf_file = vcf_sample_file(
            parameters, sample_id, ext=out_indels_ext
        )
        vcf_write_variants_features(
            out_indels_vcf_file, aggregated_indels_list
        )
    else:
        # Writing indel calls in a VCF file
        aggregated_indels_list_2 = [
            (variant, features, score_dict.copy())
            for (variant, features) in aggregated_indels_list_1
        ]
        aggregated_indels_list = get_control_samples_feature(
            parameters, aggregated_indels_list_2, control_calls
        )
        # Writing features in a temporary VCF file
        out_indels_ext = f"{FILTERED_INDELS_VCF_EXT}_aux"
        out_indels_vcf_file = vcf_sample_file(
            parameters, sample_id, ext=out_indels_ext
        )
        vcf_write_variants_features(
            out_indels_vcf_file, aggregated_indels_list
        )
        # Reading the VCF file into a DataFrame
        indels_df = vcf_import_df(out_indels_vcf_file)
        # Bug: the initial VCF files
        indels_df[ID_COL] = indels_df.apply(lambda row: '.', axis=1)
        os.remove(out_indels_vcf_file)
        indels_df.sort_values(by=[CHR_COL, POS_COL], inplace=True)
        indels_df.reset_index(drop=True, inplace=True)
        # Updating the penalty scores
        add_complexity_score(
            indels_df, manifest, COMP_KMIN, COMP_KMAX, COMP_L
        )
        add_support_score(indels_df)
        add_overlap_score(indels_df)
        add_control_score(indels_df)
        add_confidence_score(indels_df)
        indels_df[INFO_COL] = indels_df.apply(
            lambda row: reformat_info(row), axis=1
        )
        indels_df.drop(
            columns=[SCORE, COMPLEXITY, SUPPORT, OVERLAP, CONTROL],
            inplace=True
        )
        vcf_write_df(
            out_vcf_file, indels_df, sorting=VCF_SORT_POS, append=True
        )


if __name__ == "__main__":
    ARGS_RUNS = ['runs', None, 'Run ID and name csv file']
    ARGS_MANIFEST = ['-m', '--manifest','Manifest file name']

    parser = argparse.ArgumentParser(description='Indels pipeline: filters')
    parser.add_argument(ARGS_RUNS[0], type=str, help=ARGS_RUNS[2])
    parser.add_argument(ARGS_MANIFEST[0],
                        ARGS_MANIFEST[1],
                        type=str,
                        default=MANIFEST_V4_FILE,
                        help=ARGS_MANIFEST[2])
    args = parser.parse_args()

    # Loading ghe amplicons manifest and the list of non-excluded non-MSI
    # samples
    # By default MSI samples are not considered for indels calling.
    manifest_file_path = os.path.join('data', args.manifest)
    manifest = AmpliconsManifest(manifest_file_path)

    # Reading run infos
    runs_df = pd.read_csv(args.runs, header=None, names=[RUN_NAME, RUN_ID])
    for _, row in runs_df.iterrows():
        run_id = row[RUN_ID]
        run_name = row[RUN_NAME]
        print(f"{run_id}\t{run_name}")
        # Reading pipeline parameters
        parameters_file_path = os.path.join('results', run_id, f"{run_id}.yaml")
        parameters = Parameters(parameters_file_path)
        # List of samples to consider
        indels_vcf_ext_full = INDELS_VCF_EXT + '.vcf'
        sample_list = []
        for x in os.listdir(parameters.get_output_path()):
            if x.endswith(indels_vcf_ext_full):
                sample_list.append(x.replace(indels_vcf_ext_full, ''))
        # List of amplicons to consider
        amplicon_list = manifest.get_amplicons_id()
        # Filtering out MSI amplicons
        msi_amplicon_list = manifest.get_msi_amplicons_id()
        nonmsi_amplicon_list = [
            x for x in amplicon_list if x not in msi_amplicon_list
        ]

        out_vcf_file = os.path.join('results', run_id, f"{run_id}_indels.vcf")
        out_vcf = open(out_vcf_file, 'w')
        out_vcf.write(vcf_header())
        out_vcf.close()

        # Processing samples
        # Reading variants for control samples and filtering for all filters
        # but last one
        control_samples = [CONTROL_NF, CONTROL_BLANK]
        control_calls = {}
        for sample_id in sample_list:
            if parameters.check_is_control(
                sample_id, control_keys=control_samples
            ):
                print(f"\t{sample_id}")
                process_sample(
                    parameters,
                    sample_id,
                    run_id,
                    run_name,
                    nonmsi_amplicon_list,
                    manifest,
                    control_calls,
                    None,
                    control_sample=True
                )
        # Filtering variants from non-control samples
        for sample_id in sample_list:
            if not parameters.check_is_control(
                sample_id, control_keys=control_samples
            ):
                print(f"\t{sample_id}")
                process_sample(
                    parameters,
                    sample_id,
                    run_id,
                    run_name,
                    nonmsi_amplicon_list,
                    manifest,
                    control_calls,
                    out_vcf_file,
                    control_sample=False
                )
        run_snpeff(parameters, out_vcf_file)
        out_snpeff_vcf_file = out_vcf_file.replace('.vcf', '_snpeff.vcf')
        out_tsv_file = os.path.join('results', run_id, f"{run_id}_indels.tsv")
        dump_vcf_to_tsv(out_snpeff_vcf_file, out_tsv_file)
        txt_files_path = os.path.join('results', run_id)
        files_to_remove = glob.glob(f"results/{run_id}/*.txt")
        for f in files_to_remove:
            try:
                os.remove(f)
            except OSError as e:
                print("Error: %s : %s" % (f, e.strerror))
