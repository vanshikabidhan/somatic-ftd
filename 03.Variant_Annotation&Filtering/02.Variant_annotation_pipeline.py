#!/usr/bin/env python3
"""
Variant Annotation Pipeline
Create a vcf file from an h5 output, runs variant annotation using VEP, and a custom script to convert VEP CSQ to TSV
"""
import argparse
import gzip
import subprocess
import sys
from typing import Dict, List
import os
import shutil
import missionbio.mosaic as ms
import numpy as np
import pandas as pd
 
# ============================================================================
# VEP configuration and utilities
# ============================================================================

# VEP data paths configuration
VEP_CONFIG = {
    "study_path": "/results/rr/study/hg38s/study232-missionbio_TDP-C",
    "input_vcf": "/results/rr/study/hg38s/study232-missionbio_TDP-C/data/outputs/VCF_files/hg38/combined_samples_allvars_vafref10_vafhet30_vafhom90_cleaned.vcf.gz",
    "vep_output": "/results/rr/study/hg38s/study232-missionbio_TDP-C/data/outputs/annotations/hg38/combined_samples_allvars_vafref10_vafhet30_vafhom90_hgvs_annotations.vcf",
    
    # Reference data paths
    "vep_fasta": os.path.expanduser("~/.vep/homo_sapiens/115_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"),
    "loftee_dir": os.path.expanduser("~/.vep/Plugins/loftee"),
    "lof_data": os.path.expanduser("~/.vep/vep_data/loftee"),
    "alphamissense_tsv": os.path.expanduser("~/.vep/vep_data/AlphaMissense/AlphaMissense_hg38.norm.tsv.gz"),
    "cadd_snv": os.path.expanduser("~/.vep/vep_data/CADD/v1.7/whole_genome_SNVs.tsv.gz"),
    "cadd_indel": os.path.expanduser("~/.vep/vep_data/CADD/v1.7/gnomad.genomes.r4.0.indel.tsv.gz"),
    "spliceai_snv": os.path.expanduser("~/.vep/vep_data/SpliceAI/spliceai_scores.masked.snv.hg38.vcf.gz"),
    "spliceai_indel": os.path.expanduser("~/.vep/vep_data/SpliceAI/spliceai_scores.masked.indel.hg38.vcf.gz"),
}

# ============================================================================
# Create VCF from h5 output
# ============================================================================

# GRCh38 primary contigs for VCF header
HG38_CONTIGS = [
    ("chr1", 248956422), ("chr2", 242193529), ("chr3", 198295559),
    ("chr4", 190214555), ("chr5", 181538259), ("chr6", 170805979),
    ("chr7", 159345973), ("chr8", 145138636), ("chr9", 138394717),
    ("chr10", 133797422), ("chr11", 135086622), ("chr12", 133275309),
    ("chr13", 114364328), ("chr14", 107043718), ("chr15", 101991189),
    ("chr16", 90338345),  ("chr17", 83257441),  ("chr18", 80373285),
    ("chr19", 58617616),  ("chr20", 64444167),  ("chr21", 46709983),
    ("chr22", 50818468),  ("chrX", 156040895),  ("chrY", 57227415),
    ("chrM", 16569),
]

def create_vcf_from_h5(h5_path: str, output_vcf_path: str) -> None:
    """
    Create a VCF file from an h5 output (for running VEP)
    
    Args:
        h5_path: Path to the h5 file from missionbio mosaic
        output_vcf_path: Path where the VCF file will be written
    """
    print(f"Loading h5 file from: {h5_path}")
    data = ms.load(h5_path, raw=False, filter_variants=False, filter_cells=False, single=True)
    variants = data.dna.ids()

    # Prepare VCF header
    vcf_header = (
        ['##fileformat=VCFv4.2'] +
        [f'##contig=<ID={c},length={l}>' for c, l in HG38_CONTIGS] +
        ['#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO']
    )

    # Parse variants and prepare VCF lines
    vcf_lines = []
    for entry in variants:
        chrom, pos, alleles = entry.split(':')
        ref, alt = alleles.split('/')
        vcf_line = f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.'
        vcf_lines.append(vcf_line)

    # Write VCF file
    with open(output_vcf_path, 'w') as f:
        for line in vcf_header:
            f.write(line + '\n')
        for line in vcf_lines:
            f.write(line + '\n')

    print(f"VCF file created successfully at {output_vcf_path}!")
    print(f"Total variants: {len(vcf_lines)}")

# ============================================================================
# Clean VCF file
# ============================================================================

def clean_vcf(input_vcf: str, output_vcf: str) -> None:
    """
    Filter out symbolic alleles and complex variants that may cause issues with VEP.
    Uses zcat for decompression (if gzipped) and awk for filtering.
    Args:
        input_vcf: Input VCF file path (can be gzipped)
        output_vcf: Output cleaned VCF file path
    """
    print(f"Cleaning VCF file...")
    print(f"Input: {input_vcf}")
    print(f"Output: {output_vcf}")
    
    # Determine if input is gzipped
    if input_vcf.endswith('.gz'):
        # Use zcat for gzipped files
        zcat_cmd = f"zcat {input_vcf}"
    else:
        # Use cat for uncompressed files
        zcat_cmd = f"cat {input_vcf}"
    
    full_cmd = (
        f"{zcat_cmd} | "
        f"awk 'BEGIN{{FS=OFS=\"\\t\"}} /^#/ {{print; next}} "
        f"($4 ~ /[\\*\\+]/ || $5 ~ /[\\*\\+]/) {{next}} {{print}}' > {output_vcf}"
    )
    print(f"Running: {full_cmd}\n")
    
    try:
        subprocess.run(full_cmd, shell=True, check=True)
        
        # Count variants in output
        with open(output_vcf, 'r') as f:
            variant_count = sum(1 for line in f if not line.startswith('#'))
        
        print(f"VCF cleaning completed successfully!")
        print(f"Cleaned variants: {variant_count}")
    except subprocess.CalledProcessError as e:
        print(f"Error cleaning VCF: {e}")
        sys.exit(1)
    
# ============================================================================
# Run VEP 
# ============================================================================

def run_vep(input_vcf: str = None, output_vcf: str = None, config: Dict = None) -> None:
    """
    Run VEP annotation with comprehensive plugins and databases.
    
    Uses the following plugins:
    - LoF (LOFTEE)
    - AlphaMissense
    - CADD
    - SpliceAI
    
    Args:
        input_vcf: Input VCF file (default from config)
        output_vcf: Output VCF file (default from config)
        config: Custom config dict (default to VEP_CONFIG)
    """
    if config is None:
        config = VEP_CONFIG
    
    input_vcf = input_vcf or config["input_vcf"]
    output_vcf = output_vcf or config["vep_output"]
    
    vep_cmd = [
        "vep",
        "--input_file", input_vcf,
        "--output_file", output_vcf,
        "--vcf",
        "--no_stats",
        "--offline",
        "--force_overwrite",
        "--cache",
        "--dir_cache", os.path.expanduser("~/.vep"),
        "--species", "homo_sapiens",
        "--assembly", "GRCh38",
        "--fasta", config["vep_fasta"],
        "--symbol",
        "--hgvs",
        "--canonical",
        "--biotype",
        "--uniprot",
        "--mane",
        "--ccds",
        "--polyphen", "p",
        "--sift", "p",
        "--dir_plugins", config["loftee_dir"],
        "--plugin", f"LoF,loftee_path:{config['loftee_dir']},gerp_bigwig:{config['lof_data']}/gerp_conservation_scores.homo_sapiens.GRCh38.bw,human_ancestor_fa:{config['lof_data']}/human_ancestor.fa.gz,conservation_file:{config['lof_data']}/loftee.sql",
        "--plugin", f"AlphaMissense,file={config['alphamissense_tsv']},cols=all,transcript_match=0",
        "--plugin", f"CADD,snv={config['cadd_snv']},indels={config['cadd_indel']}",
        "--plugin", f"SpliceAI,snv={config['spliceai_snv']},indel={config['spliceai_indel']}",
        "--af",
        "--af_1kg",
        "--max_af",
        "--af_gnomade",
        "--af_gnomadg"
    ]
    
    print(f"Running VEP with command:")
    print(" ".join(vep_cmd))
    print()
    
    try:
        result = subprocess.run(vep_cmd, check=True)
        print(f"\nVEP annotation completed successfully!")
        print(f"Output written to: {output_vcf}")
    except subprocess.CalledProcessError as e:
        print(f"Error running VEP: {e}")
        sys.exit(1)

# ============================================================================
#  VEP CSQ field definitions and utilities
# ============================================================================

CSQ_FIELDS: List[str] = [
    "Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature","BIOTYPE",
    "EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position","Protein_position",
    "Amino_acids","Codons","Existing_variation","DISTANCE","STRAND","FLAGS",
    "SYMBOL_SOURCE","HGNC_ID","CANONICAL","MANE","MANE_SELECT","MANE_PLUS_CLINICAL",
    "CCDS","SWISSPROT","TREMBL","UNIPARC","UNIPROT_ISOFORM","SIFT","PolyPhen",
    "HGVS_OFFSET","AF","AFR_AF","AMR_AF","EAS_AF","EUR_AF","SAS_AF",
    "gnomADe_AF","gnomADe_AFR_AF","gnomADe_AMR_AF","gnomADe_ASJ_AF","gnomADe_EAS_AF",
    "gnomADe_FIN_AF","gnomADe_MID_AF","gnomADe_NFE_AF","gnomADe_REMAINING_AF","gnomADe_SAS_AF",
    "gnomADg_AF","gnomADg_AFR_AF","gnomADg_AMI_AF","gnomADg_AMR_AF","gnomADg_ASJ_AF",
    "gnomADg_EAS_AF","gnomADg_FIN_AF","gnomADg_MID_AF","gnomADg_NFE_AF","gnomADg_REMAINING_AF",
    "gnomADg_SAS_AF","MAX_AF","MAX_AF_POPS","CLIN_SIG","SOMATIC","PHENO",
    "LoF","LoF_filter","LoF_flags","LoF_info",
    "am_class","am_genome","am_pathogenicity","am_protein_variant","am_transcript_id",
    "am_uniprot_id","CADD_PHRED","CADD_RAW",
    "SpliceAI_pred_DP_AG","SpliceAI_pred_DP_AL","SpliceAI_pred_DP_DG","SpliceAI_pred_DP_DL",
    "SpliceAI_pred_DS_AG","SpliceAI_pred_DS_AL","SpliceAI_pred_DS_DG","SpliceAI_pred_DS_DL",
    "SpliceAI_pred_SYMBOL",
]

CORE_VCF_FIELDS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]

def open_gzip(path: str):
    """Open file, handling both gzip and plain text files"""
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")

def parse_info(info_str: str) -> Dict[str, str]:
    """
    Parse INFO field into dict; flag-only keys become key -> "1".
    """
    d: Dict[str, str] = {}
    if not info_str or info_str == ".":
        return d
    for item in info_str.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
        else:
            d[item] = "1"
    return d

def safe_cell(x: str) -> str:
    """Remove problematic characters from cell values"""
    return x.replace("\t", " ").replace("\n", " ").replace("\r", " ")

def explode_csq(csq_value: str) -> List[List[str]]:
    """
    Explode CSQ value into rows of annotations.
    CSQ value:
      - multiple transcript annotations separated by ','
      - fields within each annotation separated by '|'
    Returns list of lists (aligned to CSQ_FIELDS length).
    """
    if not csq_value:
        return []

    rows: List[List[str]] = []
    for ann in csq_value.split(","):        
        parts = ann.split("|")             
        if len(parts) < len(CSQ_FIELDS):
            parts += [""] * (len(CSQ_FIELDS) - len(parts))
        elif len(parts) > len(CSQ_FIELDS):
            parts = parts[:len(CSQ_FIELDS)]
        rows.append(parts)
    return rows

def convert_vcf_to_tsv(input_vcf: str, output_tsv: str, csq_key: str = "CSQ",
                       keep_info: str = "", include_empty: bool = False) -> None:
    """
    Convert VEP annotated VCF to TSV format.
    
    Args:
        input_vcf: Input VCF file path (supports .gz files, use '-' for stdin)
        output_tsv: Output TSV file path
        csq_key: INFO key containing VEP CSQ (default: CSQ)
        keep_info: Comma-separated INFO keys to include as columns (e.g. 'AC,AN,DP')
        include_empty: If CSQ is missing, still output one row with empty CSQ columns
    """
    keep_info_keys = [k.strip() for k in keep_info.split(",") if k.strip()]
    header = CORE_VCF_FIELDS + keep_info_keys + CSQ_FIELDS

    with open_gzip(input_vcf) as fin, open(output_tsv, "wt", encoding="utf-8") as fout:
        fout.write("\t".join(header) + "\n")

        for line in fin:
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                continue

            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue

            chrom, pos, vid, ref, alt, qual, flt, info_str = cols[:8]
            info = parse_info(info_str)

            base = {
                "CHROM": chrom,
                "POS": pos,
                "ID": vid,
                "REF": ref,
                "ALT": alt,
                "QUAL": qual,
                "FILTER": flt,
            }
            for k in keep_info_keys:
                base[k] = info.get(k, "")

            csq_value = info.get(csq_key, "")
            csq_rows = explode_csq(csq_value)

            if not csq_rows:
                if include_empty:
                    out = [base.get(k, "") for k in CORE_VCF_FIELDS]
                    out += [base.get(k, "") for k in keep_info_keys]
                    out += [""] * len(CSQ_FIELDS)
                    fout.write("\t".join(safe_cell(v) for v in out) + "\n")
                continue

            for parts in csq_rows:
                out = [base.get(k, "") for k in CORE_VCF_FIELDS]
                out += [base.get(k, "") for k in keep_info_keys]
                out += parts
                fout.write("\t".join(safe_cell(v) for v in out) + "\n")

    print(f"TSV conversion completed! Output written to {output_tsv}")

def main():
    ap = argparse.ArgumentParser(
        description="Variant Annotation Pipeline: Create VCF from h5, clean VCF, run VEP annotation, and convert to TSV"
    )
    subparsers = ap.add_subparsers(dest="command", help="Available commands")

    # create-vcf
    vcf_parser = subparsers.add_parser("create-vcf", help="Create VCF from h5 file")
    vcf_parser.add_argument("-i", "--input", required=True, help="Input h5 file")
    vcf_parser.add_argument("-o", "--output", required=True, help="Output VCF file")

    # clean-vcf
    clean_parser = subparsers.add_parser("clean-vcf", help="Clean VCF by removing symbolic alleles (* or +)")
    clean_parser.add_argument("-i", "--input", required=True, help="Input VCF file (can be gzipped)")
    clean_parser.add_argument("-o", "--output", required=True, help="Output cleaned VCF file")

    # run-vep
    vep_parser = subparsers.add_parser("run-vep", help="Run VEP annotation")
    vep_parser.add_argument("-i", "--input", default=None, help="Input VCF file (uses config default if not provided)")
    vep_parser.add_argument("-o", "--output", default=None, help="Output VCF file (uses config default if not provided)")

    # vcf-to-tsv
    tsv_parser = subparsers.add_parser("vcf-to-tsv", help="Convert VEP annotated VCF to TSV")
    tsv_parser.add_argument("-i", "--input", required=True, help="Input VCF/VCF.GZ (use '-' for stdin)")
    tsv_parser.add_argument("-o", "--output", required=True, help="Output TSV path")
    tsv_parser.add_argument("--csq-key", default="CSQ", help="INFO key containing VEP CSQ (default: CSQ)")
    tsv_parser.add_argument("--keep-info", default="", help="Comma-separated INFO keys to include as columns")
    tsv_parser.add_argument("--include-empty", action="store_true",
                            help="If CSQ is missing, still output one row with empty CSQ columns")

    # full-pipeline
    full_parser = subparsers.add_parser("full-pipeline", help="Run complete pipeline: h5->VCF->Clean->VEP->TSV")
    full_parser.add_argument("-h5", "--h5-input", required=True, help="Input h5 file")
    full_parser.add_argument("-vcf", "--vcf-output", required=True, help="Output VCF file")
    full_parser.add_argument("-tsv", "--tsv-output", required=True, help="Output TSV file")
    full_parser.add_argument("--skip-vep", action="store_true", help="Skip VEP annotation step")

    args = ap.parse_args()

    if args.command == "create-vcf":
        create_vcf_from_h5(args.input, args.output)
    elif args.command == "clean-vcf":
        clean_vcf(args.input, args.output)
    elif args.command == "run-vep":
        run_vep(args.input, args.output)
    elif args.command == "vcf-to-tsv":
        convert_vcf_to_tsv(args.input, args.output, args.csq_key, args.keep_info, args.include_empty)
    elif args.command == "full-pipeline":
        print("=" * 60)
        print("STEP 1: Creating VCF from h5")
        print("=" * 60)
        create_vcf_from_h5(args.h5_input, args.vcf_output)
        
        vcf_base = args.vcf_output.rsplit('.', 1)[0] if '.' in args.vcf_output else args.vcf_output
        cleaned_vcf = f"{vcf_base}_cleaned.vcf.gz"
        
        print("\n" + "=" * 60)
        print("STEP 2: Cleaning VCF (removing symbolic alleles)")
        print("=" * 60)
        # Note: The original VCF from h5 is uncompressed, so we need to handle that
        # First gzip the original VCF if it's not already
        if not args.vcf_output.endswith('.gz'):
            vcf_gz = args.vcf_output + '.gz'
            print(f"Compressing VCF: {args.vcf_output} -> {vcf_gz}")
            with open(args.vcf_output, 'rb') as f_in:
                with gzip.open(vcf_gz, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            clean_vcf(vcf_gz, cleaned_vcf)
        else:
            clean_vcf(args.vcf_output, cleaned_vcf)
        
        if not args.skip_vep:
            print("\n" + "=" * 60)
            print("STEP 3: Running VEP annotation")
            print("=" * 60)
            run_vep(cleaned_vcf, args.vcf_output)
        
        print("\n" + "=" * 60)
        print("STEP 4: Converting VCF to TSV")
        print("=" * 60)
        convert_vcf_to_tsv(args.vcf_output, args.tsv_output)
        
        print("\n" + "=" * 60)
        print("PIPELINE COMPLETE!")
        print("=" * 60)
    else:
        ap.print_help()

if __name__ == "__main__":
    main()