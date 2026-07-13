#!/bin/bash

INPUT_DIR="/results/rr/study/hg38s/study232-missionbio_TDP-C/tapestri_run/hg38/config_files_tapestri_run"
OUT_DIR="/results/rr/study/hg38s/study343-tmp_vbidhan/tapestri_v2_runs/hg38"

for cfg in "${INPUT_DIR}"/config_pool*.yaml; do
    # Skip if no matching files
    [[ -e "$cfg" ]] || { echo "No config_pool*.yaml files found"; break; }

    # Remove directory + .yaml suffix
    base=$(basename "$cfg" .yaml)
    pool=${base#config_}

    echo "Running Tapestri with config: $cfg → output-folder: $pool"

    # Run tapestri dna run
    tapestri dna run --config "$cfg" --output-folder "$OUT_DIR/$pool" --overwrite --n-cores 80
done

