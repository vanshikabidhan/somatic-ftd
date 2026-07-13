#!/bin/bash
BASE="/ori/nucl/20231016_AV224503_4606_6fYWK3D6XraZKFs"
OUT_DIR="/home/vanshika/study343-tmp_vbidhan/tapestri_v2_runs/bams"

for i in {4..25}; do
    IN_FOLDER="${BASE}/P${i}"
    OUT_FOLDER="${OUT_DIR}/P${i}"

    if [[ ! -d "$IN_FOLDER" ]]; then
        echo "Skipping P${i} (input folder not found: $IN_FOLDER)"
        continue
    fi

    echo "Processing $IN_FOLDER → $OUT_FOLDER"

    mkdir -p "$OUT_FOLDER"

    # Find R1/R2 files
    shopt -s nullglob
    R1_FILES=("$IN_FOLDER"/*_R1_001.fastq.gz)
    R2_FILES=("$IN_FOLDER"/*_R2_001.fastq.gz)
    shopt -u nullglob

    OUT_R1="${OUT_FOLDER}/P${i}_R1_merged.fastq.gz"
    OUT_R2="${OUT_FOLDER}/P${i}_R2_merged.fastq.gz"

    # Merge R1
    if (( ${#R1_FILES[@]} == 0 )); then
        echo "  No R1 files found in $IN_FOLDER → skipping R1"
    else
        echo "  Merging ${#R1_FILES[@]} R1 files → $OUT_R1"
        cat "${R1_FILES[@]}" > "$OUT_R1"
    fi

    # Merge R2
    if (( ${#R2_FILES[@]} == 0 )); then
        echo "  No R2 files found in $IN_FOLDER → skipping R2"
    else
        echo "  Merging ${#R2_FILES[@]} R2 files → $OUT_R2"
        cat "${R2_FILES[@]}" > "$OUT_R2"
    fi

    echo
done
