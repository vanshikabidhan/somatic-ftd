#!/usr/bin/env bash

## Liftover panel files from hg19 to hg38 
set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 input_hg19.bed output_hg38.bed" >&2
  exit 2
fi

IN_BED="$1"
OUT_BED="$2"

HG19_DICT="${HG19_DICT:-/home/vanshika/study232-missionbio_TDP-C/liftover/hg19.dict}"
HG38_DICT="${HG38_DICT:-/home/vanshika/study232-missionbio_TDP-C/liftover/hg38.dict}"
CHAIN="${CHAIN:-/home/vanshika/study232-missionbio_TDP-C/liftover/hg19ToHg38.over.chain.gz}"

IN_BASE="$(basename "$IN_BED")"
OUT_DIR="$(dirname "$OUT_BED")"
REJECT_BED="${OUT_DIR}/${IN_BASE}.rejected.hg38.bed"


# Check prerequisites
command -v gatk >/dev/null 2>&1 || { echo "ERROR: gatk not found in PATH" >&2; exit 1; }
[[ -s "$IN_BED" ]]   || { echo "ERROR: input BED not found/empty: $IN_BED" >&2; exit 1; }
[[ -s "$HG19_DICT" ]]|| { echo "ERROR: hg19 dict not found/empty: $HG19_DICT" >&2; exit 1; }
[[ -s "$HG38_DICT" ]]|| { echo "ERROR: hg38 dict not found/empty: $HG38_DICT" >&2; exit 1; }
[[ -s "$CHAIN" ]]    || { echo "ERROR: chain not found/empty: $CHAIN" >&2; exit 1; }


TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT

INT19="$TMPDIR/in.hg19.interval_list"
INT38="$TMPDIR/out.hg38.interval_list"
REJ_INT="$TMPDIR/rejected.interval_list"

# BED (hg19) -> interval_list (hg19)
gatk BedToIntervalList \
  -I "$IN_BED" \
  -O "$INT19" \
  -SD "$HG19_DICT"

# Liftover interval_list hg19 -> hg38
gatk LiftOverIntervalList \
  -I "$INT19" \
  -O "$INT38" \
  -SD "$HG38_DICT" \
  -CHAIN "$CHAIN" \
  -REJECT "$REJ_INT"

# Interval_list (hg38) -> BED (hg38)
grep -v '^@' "$INT38" | \
awk 'BEGIN{OFS="\t"}{
  chrom=$1; start0=$2-1; endExclusive=$3; name=$5;
  print chrom,start0,endExclusive,name
}' > "$OUT_BED"

# Export rejects 
if [[ -s "$REJ_INT" ]]; then
  grep -v '^@' "$REJ_INT" | \
  awk 'BEGIN{OFS="\t"}{
    chrom=$1; start0=$2-1; endExclusive=$3; name=$5;
    print chrom,start0,endExclusive,name
  }' > "$REJECT_BED"
else
  : > "$REJECT_BED"
fi

echo "OK: wrote hg38 BED: $OUT_BED"
echo "OK: wrote rejected BED: $REJECT_BED"



