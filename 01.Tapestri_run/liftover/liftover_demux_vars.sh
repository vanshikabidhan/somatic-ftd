#!/usr/bin/env bash

## Liftover coordinates of demultiplexing variants from hg19 to hg38 using GATK LiftoverVcf
set -euo pipefail

IN_VCF="/home/vanshika/study232-missionbio_TDP-C/data/final/demux_vars_hg19.vcf"
OUT_VCF="/home/vanshika/study232-missionbio_TDP-C/data/final/demux_vars_hg38.vcf"
REJECT_VCF="/home/vanshika/study232-missionbio_TDP-C/data/final/rejected_vars_demux_liftover.vcf"

CHAIN="/home/vanshika/study232-missionbio_TDP-C/liftover/hg19ToHg38.over.chain.gz"
REF="/home/vanshika/study232-missionbio_TDP-C/liftover/hg38.fa.gz"

## Liftover demultiplexing variants from hg19 to hg38
gatk LiftoverVcf \
  -I "$IN_VCF" \
  -O "$OUT_VCF" \
  --CHAIN "$CHAIN" \
  -R "$REF" \
  --REJECT "$REJECT_VCF"

echo "OK: $OUT_VCF"
echo "Rejected: $REJECT_VCF"
