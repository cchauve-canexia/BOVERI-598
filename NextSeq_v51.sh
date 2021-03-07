#!/usr/bin/env bash

echo "Extract results for NextSeq v51 samples"
python bin/extract_files.py data/NextSeq_v51.csv
python bin/dump_variants.py data/NextSeq_v51.csv -m CG001v5.1_Amplicon_Manifest_Panel5.1.12_20200911.tsv
python bin/check_indels.py data/NextSeq_v51.csv data/NextSeq_blacklist_03032021.tsv > trace.txt
