#!/usr/bin/env bash
# Generate config files for Tapestri runs for pools P4 to P25
for i in {4..25}; do
  P="P${i}"
  outfile="/results/rr/study/hg38s/study232-missionbio_TDP-C/tapestri_run/hg38/config_files_tapestri_run/config_pool${i}.yaml"

  cat > "${outfile}" <<EOF
flt3:
  enabled: false
gatk:
  HaplotypeCaller:
    max-reads-per-alignment-start: 100000
genome:
  path: "$HOME/tapestri/hg38/hg38.fa"
  version: hg38
inputs:
  tube1:
    r1: ["/results/rr/study/hg38s/study343-tmp_vbidhan/tapestri_v2_runs/bams/${P}/${P}_R1_merged.fastq.gz"]
    r2: ["/results/rr/study/hg38s/study343-tmp_vbidhan/tapestri_v2_runs/bams/${P}/${P}_R2_merged.fastq.gz"]
output:
  prefix: pool${i}-v2
panel:
  name: custom-panel-3971
  path: "$HOME/tapestri/panels/design-3971-hg38/"
EOF

  echo "Generated ${outfile}"
done
