#!/usr/bin/env bash

echo "Extract results for NextSeq v51 samples"
python bin/extract_files.py data/NextSeq_v51_05032021.csv
python bin/dump_variants.py data/NextSeq_v51_05032021.csv -m CG001v5.1_Amplicon_Manifest_Panel5.1.12_20200911.tsv
