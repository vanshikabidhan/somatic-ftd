# somatic-ftd

Repository for analysis of somatic variants in FTLD-TDP type C (based on the Mission Bio Tapestri platform).

## Overview

This repository contains the analysis pipeline used to call, annotate, and characterize somatic variants across a targeted gene panel from single-nucleus amplicon sequencing (Mission Bio Tapestri) of NeuN+ nuclei isolated from the superior temporal gyrus of FTLD-TDP type C patients and neuropathological controls. The pipeline comprises four main steps, corresponding to the numbered top-level directories:

1. **Run setup** — data preparation and launching the Tapestri sequencing pipeline
2. **Demultiplexing** — assign single-cell barcodes to donors from each pooled run
3. **Variant annotation & filtering** — integrate pools, annotate variants, and apply QC for selecting high-confidence variants
4. **Downstream analyses** — generate manuscript figures and supplemental analyses

## Repository Structure

```
somatic-ftd/∂
├── 01.Tapestri_run/
│   ├── data_prep/
│   │   ├── make_configs.sh
│   │   └── merge_bams.sh
│   │   └── liftover
│   │       ├── liftover_demux_vars.sh
│   │       └── liftover_panel_files.sh
│   └── run_tapestri.sh
├── 02.Demultiplexing/
│   └── demux_tapestri.py
├── 03.Variant_Annotation&Filtering/
│   ├── 01.Multisample_integration.py
│   ├── 02.Variant_annotation_pipeline.py
│   └── 03.Variant_filtering.py
├── 04.Downstream_analyses/
│   ├── Figures/
│   │   ├── Figure1.ipynb
│   │   ├── Figure2.ipynb
│   │   ├── Figure3.ipynb
│   │   └── Burden/
│   │       ├── data-prep-burden-testing.ipynb
│   │       ├── gene-model.R
│   │       └── per-domain.R
│   └── Supplemental/
│       ├── coverage-plots.ipynb
│       ├── depth-vs-vaf-plot.py
│       ├── gene-exp-scRNA-ref-MTG.R
│       ├── mosaic_loss_of_chrY.ipynb
│       └── Mutational_profiles.ipynb
├── LICENSE
└── README.md
```

## Directory & File Descriptions

### `01.Tapestri_run/` — Sequencing run setup

Scripts to prepare inputs for and launch the Tapestri DNA pipeline across sequencing pools.

- **`data_prep/make_configs.sh`** — Generates a per-pool `config_poolN.yaml` file for the `tapestri dna run` command, pointing to each pool's merged FASTQs and the hg38 reference.
- **`data_prep/merge_bams.sh`** — Merges/concatenates per-lane R1/R2 FASTQ files for each pool into a single pair of FASTQs used as input to the Tapestri run.
- **`data_prep/liftover/liftover_demux_vars.sh`** — Lifts the demultiplexing variant VCF from hg19 to hg38 using GATK `LiftoverVcf`.
- **`data_prep/liftover/liftover_panel_files.sh`** — Lifts a Tapestri panel BED file from hg19 to hg38.
- **`run_tapestri.sh`** — Iterates over all generated pool configs and runs `tapestri dna run` for each, producing per-pool `.h5` output files.

### `02.Demultiplexing/` — Sample assignment

- **`demux_tapestri.py`** — Uses Mission Bio Mosaic to demultiplex pooled single-cell data into individual samples/donors. Loads demultiplexing-variant genotypes per pool, reduces dimensionality (PCA) and clusters barcodes (KMeans) to assign each cell to a donor, and flags doublets/unassigned cells.

### `03.Variant_Annotation&Filtering/` — Integration, annotation, and QC

- **`01.Multisample_integration.py`** — Loads per-pool `.h5` files, removes doublets/unassigned cells, attaches sample metadata (e.g. group, sex, site, age at onset), splits each pool into per-sample groups, filters somatic variants, and merges all samples into a single combined `.h5`.
- **`02.Variant_annotation_pipeline.py`** — Extracts the variants from the merged `.h5` to VCF, cleans the VCF, and runs Ensembl VEP annotation (with LOFTEE, AlphaMissense, CADD, and SpliceAI plugins) before converting the annotated VCF (CSQ field) into a flat TSV. Supports running individual steps (`vcf`, `clean`, `vep`, `tsv`) or the `full` pipeline via subcommands.
- **`03.Variant_filtering.py`** — Post-integration QC and variant filtering: based on callability per sample, population frequencies, overlap with repeat regions, removal of germline variants, correction of unbalanced heterozygous calls, and masking of variants too close to amplicon boundaries.

### `04.Downstream_analyses/` — Manuscript figures and supplemental analyses

#### `Figures/`

- **`Figure1.ipynb`** — Study design and overview of all coding somatic variants: per-sample cell counts, varaint-reccurence across samples, variant frequencies.
- **`Figure2.ipynb`** — Per-codon variant distribution for TARDBP, variant frequecnies per sample, including radial plots highlighting the TARDBP p.R42H and p.E58G hotspots (Figure 2).
- **`Figure3.ipynb`** — Frequency checks for known germline TARDBP mutations and selection of variants overrepresented and exclusively present in patients (Figure 3).
- **`Burden/`** — Gene and protein-domain burden testing.
  - **`data-prep-burden-testing.ipynb`** — Prepares input matrices for gene-level and TDP-43-domain-specific burden testing of coding somatic variants.
  - **`gene-model.R`** — Gene-level burden association testing using mixed-effects models.
  - **`per-domain.R`** — TDP-43 protein-domain-specific burden association testing.

#### `Supplemental/`

- **`coverage-plots.ipynb`** — Coverage across all coding positions in STG samples post variant filtering.
- **`depth-vs-vaf-plot.py`** — Extracts heterozygous genotype calls across samples and produces a KDE marginal plot of sequencing depth vs. variant allele frequency.
- **`gene-exp-scRNA-ref-MTG.R`** — Neuronal gene-expression reference analysis using the SEA-AD snRNA-seq dataset (middle temporal gyrus), computing expression deciles for candidate genes.
- **`mosaic_loss_of_chrY.ipynb`** — Detects mosaic loss of chromosome Y (LOY) from raw read depth and Mosaic's normalization algorithm, and identifies barcodes with no chrY coverage.
- **`Mutational_profiles.ipynb`** — Mutational profiles with correction for amplicon bias

## License

Released under the [MIT License](LICENSE).