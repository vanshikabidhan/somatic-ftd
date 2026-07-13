#!/home/vanshika/mosaic/bin/python

import argparse
import os
from pathlib import Path
import missionbio.mosaic as ms
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
from matplotlib import lines as mlines
import seaborn as sns
from typing import Tuple, Union
from scipy.stats import binom
from statsmodels.stats.multitest import multipletests
import sys
print("Packages imported.")

#Logging
def setup_logging(log_path: str):
    log_file = open(log_path, 'w')
    sys.stdout = log_file
    return log_file

# ===================================
# Helpers for frequency plot
# ===================================

def sort_AP(ids: pd.Series) -> pd.Series:
    """
    Sort order of sample ID's on the plot: Controls followed by patients, and then based on numbers
    """
    k = ids.astype(str).str.extract(r'([A-Za-z]+)(\d+)_([CP])')
    order = (k[2].map({'C': 0, 'P': 1}), k[0].str.upper(), k[1].astype(int))
    order_df = pd.DataFrame({'g': order[0], 'l': order[1], 'n': order[2]})
    return ids.iloc[order_df.sort_values(['g', 'l', 'n']).index]

def plot_variant_freq_by_sample(
    df,
    id_column='ID',
    freq_col='frequency',
    stype_col='sample_type',
    jitter=0.01,
    sample_order=None,
    ax=None,
):
    """
    Scatter plot of variant frequency per sample, coloured by sample_type.
    Returns (fig, ax, samples)
    """
    d = df.copy()

    # Normalize sample_type labels
    d[stype_col] = d[stype_col].replace(
        {'P': 'FTLD-TDP type C', 'C': 'Control'}
    ).astype(str)

    # Determine sample order
    if sample_order is None:
        samples = sort_AP(pd.Series(d[id_column].unique())).tolist()
    else:
        samples = list(sample_order)

    # Palette
    pal = {
        'FTLD-TDP type C': "#F44542"   ,
        'Control': "#92C7FF"
    }
    
    created_ax = False
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6))
        created_ax = True
    else:
        fig = ax.figure

    seen_labels = set()
    for i, sample in enumerate(samples):
        g = d[d[id_column] == sample]
        if g.empty:
            continue
        x = np.random.normal(loc=i, scale=jitter, size=len(g))
        stype = g[stype_col].iloc[0]
        color = pal.get(stype, 'lightgray')
        label = stype if stype not in seen_labels else None
        if label is not None:
            seen_labels.add(stype)
        ax.scatter(x, g[freq_col].to_numpy(), color=color, s=10, label=label)

    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=90, fontsize=9)
    ax.set_xlabel('Sample ID')
    ax.set_ylabel('Variant frequency (%)')

    if seen_labels:
        ax.legend(
            loc='upper left',
            bbox_to_anchor=(1, 1),
            frameon=False,
            title=''
        )

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.tight_layout()

    return fig, ax, samples

def plot_variant_freq_by_sample_gene(
    df,
    id_column='ID',
    freq_col='frequency',
    gene_col='SYMBOL',
    target_genes=None,        
    jitter=0.01,
    sample_order=None,
    ax=None,
    other_color='lightgray'
):
    """
    Scatter plot of variant frequency per sample, coloured by gene.
    Returns (fig, ax, samples).
    """
    if target_genes is None:
        target_genes = []

    d = df.copy()
    d[gene_col] = d[gene_col].astype(str)

    # sample order
    if sample_order is None:
        samples = sort_AP(pd.Series(d[id_column].unique())).tolist()
    else:
        samples = list(sample_order)

    # palette only for target genes present in data
    present_targets = [g for g in target_genes if g in set(d[gene_col])]
    pal = sns.color_palette("tab10", n_colors=max(1, len(present_targets)))
    gene_to_color = {g: pal[i] for i, g in enumerate(present_targets)}

    # figure/axes
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6))
    else:
        fig = ax.figure

    # plot points per sample with per-point color
    for i, sample in enumerate(samples):
        g = d[d[id_column] == sample]
        if g.empty:
            continue
        x = np.random.normal(loc=i, scale=jitter, size=len(g))
        point_colors = [gene_to_color.get(gene, other_color) for gene in g[gene_col]]
        ax.scatter(x, g[freq_col].to_numpy(), s=10, c=point_colors)

    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=90, fontsize=9)
    ax.set_xlabel('Sample ID')
    ax.set_ylabel('Variant frequency (%)')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # legend
    handles = []
    labels = []
    for g in present_targets:
        handles.append(
            mlines.Line2D([0], [0], marker='o', linestyle='', markersize=6,
                          color=gene_to_color[g])
        )
        labels.append(g)

    if (d[gene_col].isin(present_targets) == False).any():
        handles.append(
            mlines.Line2D([0], [0], marker='o', linestyle='', markersize=6,
                          color=other_color)
        )
        labels.append('Other')

    if handles:
        ax.legend(
            handles, labels,
            loc='upper left',
            bbox_to_anchor=(1, 1),
            frameon=False,
            title='Gene'
        )

    fig.tight_layout()
    return fig, ax, samples

def plot_variant_freq_by_sample_callability(
    df,
    id_column='ID',
    freq_col='frequency',
    callability_col='callability',
    jitter=0.01,
    sample_order=None,
    ax=None,
    other_color='lightgray',   # used for NaN / out-of-range
    step=10,                   # 10-point bins: 0–10, 11–20, ...
    cmap_name="viridis"      
):
    """
    Scatter plot of variant frequency per sample, coloured by callability bins.
    Callability expected in [0, 100]. Bins: 0–10, 11–20, ..., 91–100 (by default).

    Returns (fig, ax, samples).
    """
    d = df.copy()
    d[callability_col] = pd.to_numeric(d[callability_col], errors="coerce")
    # sample order
    if sample_order is None:
        samples = sort_AP(pd.Series(d[id_column].unique())).tolist()
    else:
        samples = list(sample_order)

    # Build bins: 0-10, 11-20, ..., 91-100 (inclusive)
    edges = np.arange(0, 100 + step, step)
    def _bin_label(lo, hi):
        if lo == 0:
            return f"{int(lo)}–{int(hi)}"
        return f"{int(lo)+1}–{int(hi)}"

    labels = [_bin_label(edges[i], edges[i+1]) for i in range(len(edges) - 1)]

    d["_call_bin"] = pd.cut(
        d[callability_col],
        bins=edges,
        labels=labels,
        include_lowest=True,
        right=True
    ).astype("object")  # allows NaN to stay NaN


    cmap = plt.get_cmap(cmap_name)
    pal = [cmap(i) for i in np.linspace(0.1, 0.9, len(labels))]
    bin_to_color = dict(zip(labels, pal))

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6))
    else:
        fig = ax.figure

    # plot points per sample with per-point color
    for i, sample in enumerate(samples):
        g = d[d[id_column] == sample]
        if g.empty:
            continue
        x = np.random.normal(loc=i, scale=jitter, size=len(g))
        point_colors = [
            bin_to_color.get(b, other_color) if b is not None else other_color
            for b in g["_call_bin"].to_numpy()
        ]
        ax.scatter(x, g[freq_col].to_numpy(), s=10, c=point_colors)

    ax.set_xticks(range(len(samples)))
    ax.set_xticklabels(samples, rotation=90, fontsize=9)
    ax.set_xlabel('Sample ID')
    ax.set_ylabel('Variant frequency (%)')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # legend (all bins that are present)
    present_bins = [lab for lab in labels if (d["_call_bin"] == lab).any()]
    handles, leg_labels = [], []

    for lab in present_bins:
        handles.append(
            mlines.Line2D([0], [0], marker='o', linestyle='', markersize=6,
                          color=bin_to_color[lab])
        )
        leg_labels.append(lab)

    if d["_call_bin"].isna().any():
        handles.append(
            mlines.Line2D([0], [0], marker='o', linestyle='', markersize=6,
                          color=other_color)
        )
        leg_labels.append('Missing/Out of range')

    if handles:
        ax.legend(
            handles, leg_labels,
            loc='upper left',
            bbox_to_anchor=(1, 1),
            frameon=False,
            title=f'Callability ({step}-pt bins)'
        )

    fig.tight_layout()
    return fig, ax, samples

# ====================================
# Compute variant frequency per sample
# ====================================
def compute_variant_frequencies(
    NGT_df: pd.DataFrame,
    n_meta_cols: int = 3,
) -> pd.DataFrame:
    """
    Compute per-sample variant frequency in long format. (n_rows: n_variants x n_samples )

    Parameters
    ----------
    NGT_df : pd.DataFrame
        Genotype dataframe with:
        - Variants as columns 
        - Last `n_meta_cols` columns: at least 'sample_id', 'sample_type', 'ID'
    n_meta_cols : int, default=3
        Number of non-variant columns at the end.

    Returns
    -------
    out : pd.DataFrame
        Long-format table:
        ['sample_id', 'variant_id', 'frequency', 'count',
         'n_cells_callable', 'ID', 'sample_type',
         'n_cells_total', 'callability']
    """
    variant_cols = list(NGT_df.columns[:-n_meta_cols])

    NGT_FILTERED_final = NGT_df.copy()
    sample_ids = NGT_FILTERED_final['sample_id'].unique()

    total_vars = pd.DataFrame(index=sample_ids, columns=variant_cols)
    mut_count = pd.DataFrame(index=sample_ids, columns=variant_cols)
    n_cells_genotyped = pd.DataFrame(index=sample_ids, columns=variant_cols)
    other_variables = pd.DataFrame(
        index=sample_ids,
        columns=['ID', 'sample_type', 'n_cells_total', 'max_freq', 'min_freq']
    )

    for sample in sample_ids:
        sub = NGT_FILTERED_final[NGT_FILTERED_final['sample_id'] == sample]

        sample_type = sub['sample_type'].iat[0]
        ID = sub['ID'].iat[0]
        n_cells = len(sub)

        genotyped_cells = sub[variant_cols].isin([0, 1, 2]).sum()
        counts = sub[variant_cols].isin([1, 2]).sum()
        freqs = counts / genotyped_cells * 100

        total_vars.loc[sample, variant_cols] = freqs
        mut_count.loc[sample, variant_cols] = counts
        n_cells_genotyped.loc[sample, variant_cols] = genotyped_cells

        other_variables.loc[sample, 'sample_type'] = sample_type
        other_variables.loc[sample, 'ID'] = ID
        other_variables.loc[sample, 'n_cells_total'] = n_cells
        other_variables.loc[sample, 'max_freq'] = freqs.max()
        other_variables.loc[sample, 'min_freq'] = freqs.min()

    # Long version of frequencies
    out = (
        total_vars.apply(pd.to_numeric, errors='coerce')
        .stack()
        .rename('frequency')
        .reset_index()
        .rename(columns={'level_0': 'sample_id',
                         'level_1': 'variant_id'})
    )

    # Long counts
    counts_long = (
        mut_count.apply(pd.to_numeric, errors='coerce')
        .stack()
        .rename('count')
        .reset_index()
        .rename(columns={'level_0': 'sample_id',
                         'level_1': 'variant_id'})
    )

    # Long callable cells
    cells_long = (
        n_cells_genotyped.apply(pd.to_numeric, errors='coerce')
        .stack()
        .rename('n_cells_callable')
        .reset_index()
        .rename(columns={'level_0': 'sample_id',
                         'level_1': 'variant_id'})
    )

    # Merge helper columns
    out = (
        out.merge(counts_long, on=['sample_id', 'variant_id'], how='left')
           .merge(cells_long, on=['sample_id', 'variant_id'], how='left')
    )

    ov = other_variables.rename_axis('sample_id').reset_index()

    # Add sample-level metadata
    out = out.merge(
        ov[['sample_id', 'ID', 'sample_type', 'n_cells_total']],
        on='sample_id',
        how='left'
    )

    # Callability
    out['callability'] = out['n_cells_callable'] / out['n_cells_total'] * 100

    return out

def tag_repeats(variants: pd.DataFrame, repeats: pd.DataFrame, bp=5) -> pd.DataFrame:
    """
    Tag variants with whether they fall in/near repeat regions.
    """
    def to_chr(series: pd.Series) -> pd.Series:
        s = series.astype(str).str.strip()
        s = s.str.replace(r'^(chr)?', '', regex=True)  # remove any leading 'chr'
        return 'chr' + s

    v = variants.rename(columns={"position": "pos"}).copy()
    r = repeats.rename(columns={"chrom": "chr"}).copy()

    v["chr"] = to_chr(v["chr"])
    r["chr"] = to_chr(r["chr"])

    v["pos"] = pd.to_numeric(v["pos"], errors="coerce")
    r["start"] = pd.to_numeric(r["start"], errors="coerce")
    r["end"] = pd.to_numeric(r["end"], errors="coerce")

    v = v.dropna(subset=["pos"])
    r = r.dropna(subset=["start", "end"])

    r["start5"] = (r["start"] - bp).clip(lower=1)
    r["end5"] = r["end"] + bp

    v["in_repeat"] = False
    v["near_repeat_5bp"] = False
    for col in ["repClass", "repFamily", "repName"]:
        v[col] = pd.NA

    for cchr in v["chr"].unique():
        v_mask = v["chr"] == cchr
        r_chr = r[r["chr"] == cchr].sort_values("start")
        if r_chr.empty:
            continue

        v_sorted = v.loc[v_mask].sort_values("pos")

        # inside repeat
        m_in = pd.merge_asof(
            v_sorted,
            r_chr[["start", "end", "repClass", "repFamily", "repName"]],
            left_on="pos",
            right_on="start",
            direction="backward",
            allow_exact_matches=True,
            suffixes=("", "_rep")
        ).set_index(v_sorted.index).sort_index()

        in_mask = (m_in["pos"] <= m_in["end"]).fillna(False)
        v.loc[in_mask.index, "in_repeat"] = in_mask.values
        for col in ["repClass", "repFamily", "repName"]:
            v.loc[in_mask.index, col] = m_in[f"{col}_rep"].where(in_mask).values

        # near boundary (±bp), but not inside
        r5 = r_chr[["start5", "end5"]].rename(columns={"start5": "start", "end5": "end"})
        m_near = pd.merge_asof(
            v_sorted,
            r5[["start", "end"]],
            left_on="pos",
            right_on="start",
            direction="backward",
            allow_exact_matches=True
        ).set_index(v_sorted.index).sort_index()

        near_any = (m_near["pos"] <= m_near["end"]).fillna(False)
        v.loc[near_any.index, "near_repeat_5bp"] = (
            near_any & ~v.loc[near_any.index, "in_repeat"]
        ).values

    return v

# ========================================================
# Helpers for binomial thresholding for false positives
# ========================================================

def _materialize_p0_matrix(p0: Union[float, np.ndarray, pd.DataFrame],
                           like_df: pd.DataFrame,
                           shape: Tuple[int, int]) -> np.ndarray:
    if isinstance(p0, (int, float)):
        return np.full(shape, float(p0), dtype=float)
    if isinstance(p0, pd.DataFrame):
        return p0.reindex_like(like_df).to_numpy(dtype=float, copy=False)
    p_mat = np.asarray(p0, dtype=float)
    if p_mat.shape != shape:
        raise ValueError(f"p0 array shape {p_mat.shape} must match {shape}")
    return p_mat

def pvalues_one_sided_less_vec(n: np.ndarray,
                               alt: np.ndarray,
                               het_mask: np.ndarray,
                               p0_mat: np.ndarray) -> np.ndarray:
    """p = P(K <= k | Binom(n, p0)) for het entries; NaN elsewhere."""
    pvals = np.full(n.shape, np.nan, dtype=float)
    rows, cols = np.where(het_mask)
    if rows.size == 0:
        return pvals

    N = n[rows, cols]
    K = alt[rows, cols]
    P = p0_mat[rows, cols]

    valid = (N > 0) & np.isfinite(P)
    if np.any(valid):
        Pv = np.clip(P[valid], 1e-12, 1.0 - 1e-12)
        pvals_sub = binom.cdf(K[valid], N[valid], Pv)
        pvals[rows[valid], cols[valid]] = pvals_sub
    return pvals

def convert_unbalanced_hets(
    NGT: pd.DataFrame,
    DP: pd.DataFrame,
    AF_percent: pd.DataFrame,
    *,
    alpha: float = 0.05,                 # BH-FDR threshold per cell
    p0: Union[float, np.ndarray, pd.DataFrame] = 0.5,
    min_dp: int = 10,                     # ignore hets with DP < min_dp
    apply_only_below_50: bool = True    
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Flip heterozygous calls (NGT==1) to REF (0) if:
      - one-sided binomial p-value (less) BH-FDR q < alpha
    with ALT reconstructed as round((AF%/100) * DP).
    AF_percent is assumed to be in [0,100] (or [1,100]) for sure.

    Returns: (NGT_out, flipped_report, kept_report)
    """
    DP_view = DP.reindex_like(NGT)
    AF_view = AF_percent.reindex_like(NGT)

    NGT_out  = NGT.copy()
    het_mask = (NGT_out.to_numpy(copy=False) == 1)

    # total depth
    dp = np.nan_to_num(DP_view.to_numpy(dtype=float, copy=False), nan=0.0)
    n = np.rint(dp).astype(np.int32, copy=False)
    n[n < 0] = 0

    # Convert AF% to fraction
    afp = np.nan_to_num(AF_view.to_numpy(dtype=float, copy=False), nan=0.0)
    af_frac = np.clip(afp / 100.0, 0.0, 1.0)

    if apply_only_below_50:
        het_mask = het_mask & (af_frac < 0.5)

    alt = np.rint(af_frac * dp).astype(np.int32, copy=False)
    np.clip(alt, 0, None, out=alt)
    np.minimum(alt, n, out=alt)

    if min_dp > 1:
        het_mask = het_mask & (n >= int(min_dp))

    p0_mat = _materialize_p0_matrix(p0, like_df=NGT, shape=n.shape)
    # one-sided p-values
    p_one = pvalues_one_sided_less_vec(n, alt, het_mask, p0_mat)

    # BH-FDR per cell (row-wise) across hets
    p_adj = np.full_like(p_one, np.nan, dtype=float)
    H, W = n.shape
    for i in range(H):
        row_mask = het_mask[i, :]
        if not np.any(row_mask):
            continue
        row_p = p_one[i, row_mask]
        valid = ~np.isnan(row_p)
        if not np.any(valid):
            continue
        _, qvals, *_ = multipletests(row_p[valid], alpha=alpha, method="fdr_bh")
        idx_cols = np.flatnonzero(row_mask)[valid]
        p_adj[i, idx_cols] = qvals

    flip_mask = het_mask & (p_adj < alpha)
    keep_mask = (NGT_out.to_numpy(copy=False) == 1) & ~flip_mask  # keep all other hets

    out_arr = NGT_out.to_numpy(copy=False)
    out_arr[flip_mask] = 0

    def _build_report(mask: np.ndarray, decision: str) -> pd.DataFrame:
        r, c = np.where(mask)
        if r.size == 0:
            return pd.DataFrame(columns=[
                "cell","variant","total_depth","af_percent","alt_reads_est",
                "p0_expected","p_one_sided_less","q_fdr_cell","decision"
            ])
        return pd.DataFrame({
            "cell":             NGT.index.take(r),
            "variant":          NGT.columns.take(c),
            "total_depth":      n[r, c],
            "af_percent":       afp[r, c],
            "alt_reads_est":    alt[r, c],
            "p0_expected":      p0_mat[r, c],
            "p_one_sided_less": p_one[r, c],
            "q_fdr_cell":       p_adj[r, c],
            "decision":         decision
        }).sort_values(["variant", "cell"], ignore_index=True)

    flipped_report = _build_report(flip_mask, "flipped_to_ref")
    kept_report    = _build_report(keep_mask,  "kept_as_het")

    return NGT_out, flipped_report, kept_report

def extract_pos(variant_id):
    return int(variant_id.split(":")[1])

def mark_close_variants_as_missing(geno, variant_ids, amplicon_map, max_distance=25):
    """
    variant_ids: list of variant IDs
    amplicon_map: pd.Series with index=variant_id, value=amplicon_id
    """
    n_cells, n_variants = geno.shape
    positions = np.array([extract_pos(v) for v in variant_ids])
    # Map each variant to an amplicon
    amplicon_ids = amplicon_map.reindex(variant_ids).values
    
    # Precompute amplicon → list of variant indices
    amplicon_to_variants = {}
    for i, amp in enumerate(amplicon_ids):
        amplicon_to_variants.setdefault(amp, []).append(i)
    
    # Modified genotype output
    geno_out = geno.copy()

    # Process each cell
    cell_counter = 0
    variants_counter = 0
    for s in range(n_cells):
        # genotypes for this cell
        row = geno[s, :]
        
        for amp, variant_idx_list in amplicon_to_variants.items():
            idx = np.array(variant_idx_list)

            het_hom = np.where(np.isin(row[idx], [1,2]))[0]
            if len(het_hom) < 2:
                continue

            # Check pairwise distances
            abs_positions = positions[idx][het_hom]
            diff_matrix = np.abs(abs_positions[:, None] - abs_positions[None, :])
            
            # If any pair < max_distance
            if np.any((diff_matrix < max_distance) & (diff_matrix > 0)):
                # Mark *all* those GTs as 3
                geno_out[s, idx[het_hom]] = 3
                cell_counter+=1
                variants_counter += len(het_hom)
    print(f"Total cells affected: {cell_counter}, total variants marked as missing: {variants_counter}")
    return geno_out

# =========================
# Main pipeline
# =========================

def main():
    parser = argparse.ArgumentParser(
        description="Variant filtering and quality control post integration of pools."
    )
    parser.add_argument(
        "--h5",
        default="/home/AD/vbidhan/study232-missionbio_TDP-C/data/outputs/integrated_h5_files/hg38/combined_samples_allvars_vafref10_vafhet30_vafhom90.h5",
        help="Path to integrated H5 file."
    )
    parser.add_argument(
        "--amplicons",
        default="/home/AD/vbidhan/study232-missionbio_TDP-C/design/hg38/amplicons-with-gene-annotation.tsv",
        help="Amplicon design file with gene annotation."
    )
    parser.add_argument(
        "--repeats",
        default="/home/AD/vbidhan/study232-missionbio_TDP-C/data/repeats_in_overlapping_regions_hg38.tsv",
        help="RepeatMasker overlaps file."
    )
    parser.add_argument(
        "--annotations",
        default="/home/AD/vbidhan/study232-missionbio_TDP-C/data/outputs/annotations/hg38/combined_samples_allvars_vafref10_vafhet30_vafhom90_hgvs_annotations.tsv", 
        help="Variant annotation file (VEP) " \
        "Includes output from: LOFTEE, gnomAD ferquency, SpliceAI and AlphaMissense scores."
    )
    parser.add_argument(
        "--outdir",
        default="/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/",
        help="Base output directory."
    )

    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    parquet_dir = outdir / "parquets"
    parquet_dir.mkdir(parents=True, exist_ok=True)

    log_path = outdir/ 'QC_post_integration.log'
    # --- Logging ---
    log_file = setup_logging(log_path)

    # -------------------------------------------------------
    # Load mosaic data (h5)
    # -------------------------------------------------------
    print("Loading H5 data...")
    data = ms.load(args.h5, raw=False, filter_variants=False, filter_cells=False, single=True)

    NGT_FILTERED = data.dna.get_attribute('NGT_FILTERED')
    AF_FILTERED = data.dna.get_attribute('AF_FILTERED')
    DP = data.dna.get_attribute('DP')
    GQ = data.dna.get_attribute('GQ')

    # Change Variant Identifier -> Amplicon ID
    NGT_FILTERED.columns = data.dna.col_attrs['amplicon']
    DP.columns = data.dna.col_attrs['amplicon']
    AF_FILTERED.columns = data.dna.col_attrs['amplicon']
    GQ.columns = data.dna.col_attrs['amplicon']

    # -------------------------------------------------------
    # Filter poor-performing amplicons
    # -------------------------------------------------------
    print("Filtering amplicons with good coverage across >80 percent barcodes...")
    DP_collapsed = DP.T.groupby(DP.columns).sum(min_count=1).T
    poor_amplicons = pd.DataFrame((DP_collapsed < 10).sum())
    poor_amplicons.columns = ['cells_with_avg_depth_below10']

    amplicons = pd.read_csv(
        args.amplicons,
        sep='\t',
        header=None
    )
    amplicons.columns = ['chr', 'start_pos', 'end_pos', 'amplicon_id',
                         'gene_id', 'gene_name']

    allosomes_amplicons = [
        "AMPL293806", "AMPL293807", "AMPL293808",
        "AMPL293809", "AMPL293812", "AMPL293803",
        "AMPL293804", "AMPL293805", 
        "AMPL349615", "AMPL349616", "AMPL349617", "AMPL349618"]
    amplicons.loc[amplicons['amplicon_id'].isin(allosomes_amplicons), 'gene_name'] = 'chrY'

    gene_map = (
        amplicons.assign(id=amplicons['amplicon_id'].astype(str))
        .drop_duplicates('amplicon_id')
        .set_index('amplicon_id')['gene_name']
    )

    poor_amplicons.index = poor_amplicons.index.astype(str)
    poor_amplicons['gene_name'] = poor_amplicons.index.map(gene_map)

    # Export coverage summary
    poor_amplicons.to_csv(outdir / "amplicons_depth_summary.tsv", sep="\t")

    amplicons_to_be_dropped = poor_amplicons[
        poor_amplicons['cells_with_avg_depth_below10'] >= (0.20 * data.dna.shape[0])
    ].index.tolist()

    amplicons_to_be_dropped = list(set(amplicons_to_be_dropped +  allosomes_amplicons ))
    print(f"Amplicons removed: {len(amplicons_to_be_dropped)}")
    bad_performing_aplicons = [amp for amp in amplicons_to_be_dropped if amp not in allosomes_amplicons]
    print("Bas performing amplicons removed:", bad_performing_aplicons  )

    idxs = [i for i, c in enumerate(DP.columns) if c in amplicons_to_be_dropped]
    variants_to_drop = sorted(set(i for i in idxs if 0 <= i < NGT_FILTERED.shape[1]))
    
    print("Variants dropped due to removal of sex chromsomes:",  len([i for i, c in enumerate(DP.columns) if c in allosomes_amplicons]))
    print("Variants dropped due to removal of bas performing amplicons:",  len([i for i, c in enumerate(DP.columns) if c in bad_performing_aplicons]))
    print("Total variants dropped:", len(variants_to_drop))
    NGT_FILTERED_new = NGT_FILTERED.drop(columns=NGT_FILTERED.columns[variants_to_drop])
    DP_new = DP.drop(columns=DP.columns[variants_to_drop])
    AF_FILTERED_new = AF_FILTERED.drop(columns=AF_FILTERED.columns[variants_to_drop])
    GQ_new = GQ.drop(columns=GQ.columns[variants_to_drop])

    # -------------------------------------------------------
    # Filter cells with low average coverage 
    # -------------------------------------------------------
    print("Filtering cells with good coverage across >80% amplicons")
    DP_collapsed_cleaned = DP_collapsed.drop(columns=amplicons_to_be_dropped)
    cells_avgcoverage_below10 = pd.DataFrame(
        DP_collapsed_cleaned.select_dtypes(include='number').lt(10).sum(axis=1)
    )
    cells_avgcoverage_below10.columns = ['n_amplicons_avgcoverage_below10']
    cells_avgcoverage_below10['sample_id'] = data.dna.get_attribute('sample_id')

    cells_basedonamplicons = cells_avgcoverage_below10[
        cells_avgcoverage_below10['n_amplicons_avgcoverage_below10']
        >= (0.20 * DP_collapsed_cleaned.shape[1])
    ].index

    print(f"Cells removed (low coverage >=20% amplicons): {len(cells_basedonamplicons)}")

    # Export summary - average coverage across each barcode
    cells_avgcoverage_below10.to_csv(outdir / "cells_low_coverage_amplicons_summary.tsv", sep="\t")

    # Remove cells
    NGT_FILTERED_new = NGT_FILTERED_new.drop(index=cells_basedonamplicons)
    DP_new = DP_new.drop(index=cells_basedonamplicons)
    AF_FILTERED_new = AF_FILTERED_new.drop(index=cells_basedonamplicons)
    GQ_new = GQ_new.drop(index=cells_basedonamplicons)

    # -------------------------------------------------------
    # Map back to variant IDs
    # -------------------------------------------------------
    print("Mapping back column names to variant IDs...")
    variant_ids = data.dna.get_attribute('id')
    ids_toberemoved = variant_ids.columns[variants_to_drop]
    variant_ids = variant_ids.drop(columns=ids_toberemoved)

    NGT_FILTERED_new.columns = variant_ids.columns
    DP_new.columns = variant_ids.columns
    AF_FILTERED_new.columns = variant_ids.columns
    GQ_new.columns = variant_ids.columns

    # ------------------------------------------------------- 
    # Remove variant IDs with a spanning deletion allele (marked with asterisk)
    # -------------------------------------------------------
    masked_vars = [col for col in NGT_FILTERED_new.columns if re.search(r'\*', col)]
    print("Variant IDs removed due to spanning deletion:", len(masked_vars))

    NGT_FILTERED_new.drop(masked_vars, axis=1, inplace=True)
    DP_new.drop(masked_vars, axis=1, inplace=True)
    AF_FILTERED_new.drop(masked_vars, axis=1, inplace=True)
    GQ_new.drop(masked_vars, axis=1, inplace=True)

    # -------------------------------------------------------
    # Attach metadata
    # -------------------------------------------------------
    print("Building metadata table...")
    NGT_FILTERED_new['ID'] = data.dna.get_attribute('ID').drop(index=cells_basedonamplicons)
    NGT_FILTERED_new['sample_id'] = data.dna.get_attribute('sample_id').drop(index=cells_basedonamplicons)
    NGT_FILTERED_new['sample_type'] = data.dna.get_attribute('sample_type').drop(index=cells_basedonamplicons)

    metadata_filtered = pd.DataFrame()
    metadata_filtered['barcode'] = data.dna.get_attribute('barcode').drop(index=cells_basedonamplicons)
    metadata_filtered['ID'] = data.dna.get_attribute('ID').drop(index=cells_basedonamplicons)
    metadata_filtered['age_at_death'] = data.dna.get_attribute('age_at_death').drop(index=cells_basedonamplicons)
    metadata_filtered['age_at_onset'] = data.dna.get_attribute('age_at_onset').drop(index=cells_basedonamplicons)
    metadata_filtered['sample_id'] = data.dna.get_attribute('sample_id').drop(index=cells_basedonamplicons)
    metadata_filtered['sex'] = data.dna.get_attribute('sex').drop(index=cells_basedonamplicons)
    metadata_filtered['sample_type'] = data.dna.get_attribute('sample_type').drop(index=cells_basedonamplicons)
    metadata_filtered['hemisphere'] = data.dna.get_attribute('hemisphere_ab').drop(index=cells_basedonamplicons)
    metadata_filtered['PMI'] = data.dna.get_attribute('PMI').drop(index=cells_basedonamplicons)
    metadata_filtered['doublet_score'] = data.dna.get_attribute('doublet_score').drop(index=cells_basedonamplicons)
    metadata_filtered['ethnicity'] = data.dna.get_attribute('ethnicity').drop(index=cells_basedonamplicons)
    metadata_filtered['sample_tissue'] = data.dna.get_attribute('sample_tissue').drop(index=cells_basedonamplicons)
    metadata_filtered['site'] = data.dna.get_attribute('site').drop(index=cells_basedonamplicons)
    metadata_filtered['pool'] = data.dna.get_attribute('pool').drop(index=cells_basedonamplicons)
    metadata_filtered['race'] = data.dna.get_attribute('race').drop(index=cells_basedonamplicons)
    metadata_filtered.set_index('barcode', inplace=True)

    # -------------------------------------------------------
    # Plot before downstream filtering
    # -------------------------------------------------------
    print("Plotting pre-filter variant frequencies...")
    out = compute_variant_frequencies(NGT_FILTERED_new)
    fig0, ax0, sample_order = plot_variant_freq_by_sample(out)
    fig0.savefig(figdir / "fig00_variant_freq_by_sample_pre_filters.png",  transparent=True, dpi=600)
    plt.close(fig0)

    # -------------------------------------------------------
    # Annotate variants  - repititive region 
    # -------------------------------------------------------
    print("Annotating variants with repeat information...")
    variant_cols = list(NGT_FILTERED_new.columns[:-3])
    variants_df = pd.DataFrame({'variant_id': variant_cols})
    parts = variants_df['variant_id'].astype(str).str.split(':', n=2, expand=True)
    variants_df[['chr', 'position', 'var']] = parts
    variants_df['position'] = pd.to_numeric(variants_df['position'], errors='coerce')

    # Load repeats
    rep = pd.read_csv(
        args.repeats,
        sep="\t",
        header=None,
        usecols=[0, 1, 2, 3, 5, 11, 12],
        names=["chrom", "start", "end", "repName", "strand", "repClass", "repFamily"]
    )
    rep["start"] = pd.to_numeric(rep["start"], errors="coerce").astype("Int64")
    rep["end"] = pd.to_numeric(rep["end"], errors="coerce").astype("Int64")
    rep = rep.dropna(subset=["start", "end"]).astype({"start": "int64", "end": "int64"})

    repeat_variant_df = tag_repeats(variants_df, rep, bp=5)

    # -------------------------------------------------------
    # Map variants to genes (annotations)
    # -------------------------------------------------------
    print("Loading variant annotations...")
    res = pd.read_csv(args.annotations, sep="\t")
    res['variant_id'] = res.apply(
            lambda r: f"chr{r['CHROM']}:{r['POS']}:{r['REF']}/{r['ALT']}",
            axis=1
        )
    gene_of_interest = ["TARDBP", "OPTN", "TBK1", "GRN", "UNC13A","TET2", "TMEM106B"]
    annotations = res[res["SYMBOL"].isin(gene_of_interest)]
    annotations = annotations[annotations['CANONICAL'] == "YES"]
    annotations.index = annotations['variant_id']

    variants_removed = annotations.index[~annotations.index.isin(NGT_FILTERED_new.columns[:-3])]
    variants_untargeted_regions = NGT_FILTERED_new.iloc[:, :-3].columns[
        ~NGT_FILTERED_new.iloc[:, :-3].columns.isin(annotations.index)
    ]

    NGT_FILTERED_new.drop(variants_untargeted_regions, axis=1, inplace=True)
    DP_new.drop(variants_untargeted_regions, axis=1, inplace=True)
    AF_FILTERED_new.drop(variants_untargeted_regions, axis=1, inplace=True)
    GQ_new.drop(variants_untargeted_regions, axis=1, inplace=True)
    print(len(variants_untargeted_regions), "variants from multiplexing amplicons and distant genes removed")

    annotations.drop(index=variants_removed, inplace=True)

    map_df = annotations[
        ['variant_id', 'Consequence', 'IMPACT', 'SYMBOL',
            'AF', 'gnomADe_AF', 'gnomADg_AF']
    ].set_index('variant_id')
    repeat_variant_df.set_index('variant_id', inplace=True)
    repeat_variant_df = repeat_variant_df.join(map_df, on='variant_id', how='left')
    repeat_variant_df = repeat_variant_df.loc[(NGT_FILTERED_new.columns[:-3])]
    
    # -------------------------------------------------------
    # Filter variants: common gnomAD + repeats
    # -------------------------------------------------------
    print("Filtering common and repetitive variants...")
    variants_inrepeats = repeat_variant_df[repeat_variant_df['in_repeat']].index.tolist()
    print("Variants in repeat regions", len(variants_inrepeats))

    annotations.drop(variants_inrepeats, inplace=True)
    gnomAD_vars = list(set(
        annotations[annotations['gnomADe_AF'] >= 0.01]['variant_id'].tolist()
        + annotations[annotations['gnomADg_AF'] >= 0.01]['variant_id'].tolist()
        + annotations[annotations['AF'] >= 0.01]['variant_id'].tolist()
    ))

    print("GnomAD vars:", len(gnomAD_vars))
    vars_remove_gnomad_repeats = list(set(gnomAD_vars + variants_inrepeats))
    print(f"Total variants removed (gnomAD>1% or in repeats): {len(vars_remove_gnomad_repeats)}")

    NGT_FILTERED_new.drop(vars_remove_gnomad_repeats, axis=1, inplace=True)
    DP_new.drop(vars_remove_gnomad_repeats, axis=1, inplace=True)
    AF_FILTERED_new.drop(vars_remove_gnomad_repeats, axis=1, inplace=True)
    GQ_new.drop(vars_remove_gnomad_repeats, axis=1, inplace=True)
    annotations.drop(gnomAD_vars, inplace=True)
    
    # -------------------------------------------------------
    # Recompute frequencies post-filtering
    # -------------------------------------------------------
    print("Recomputing frequencies after variant filters...")
    out = compute_variant_frequencies(NGT_FILTERED_new)
    out_filtered = out[out['frequency'] > 0].copy()

    parts = out_filtered['variant_id'].astype(str).str.split(':', n=2, expand=True)
    out_filtered[['chr', 'position', 'var']] = parts
    out_filtered['position'] = pd.to_numeric(out_filtered['position'], errors='coerce')

    gene_lookup = annotations['SYMBOL']
    out_filtered['SYMBOL'] = out_filtered['variant_id'].map(gene_lookup)

    # Plot after removing repeats + common gnomAD variants
    fig1, ax1, _ = plot_variant_freq_by_sample(out_filtered, sample_order=sample_order)
    fig1.savefig(figdir / "fig1_variant_freq_by_sample_after_repeats_gnomad.png",  transparent=True, dpi=600)
    plt.close(fig1)

    # Callability plot 
    fig2, ax2, _ = plot_variant_freq_by_sample_callability(out_filtered, sample_order=sample_order)
    fig2.savefig(figdir / "fig2_variant_freq_by_callability.png",  transparent=True, dpi=600)
    plt.close(fig2)
    
    # -------------------------------------------------------
    # Filter on callability >=70
    # -------------------------------------------------------
    print("Filtering on callability >=70%...")
    out_filter1 = out_filtered[out_filtered['callability'] >= 70].copy()
    fig3, ax3, _ = plot_variant_freq_by_sample(out_filter1, sample_order=sample_order)
    fig3.savefig(figdir / "fig3_variant_freq_by_sample_callability70.png",  transparent=True, dpi=600)
    plt.close(fig3)

    # Per gene plot 
    fig4, ax4, _ = plot_variant_freq_by_sample_gene(
        out_filter1,
        id_column='ID',
        freq_col='frequency',
        gene_col='SYMBOL',
        target_genes=gene_of_interest,
        sample_order=sample_order
    )
    fig4.savefig(figdir / "fig4_variant_freq_by_sample_gene_callability70.png",  transparent=True, dpi=600)
    plt.close(fig4)

    # -------------------------------------------------------
    # Remove putative germline variants (>=65% freq)
    # -------------------------------------------------------
    print("Tagging and removing putative germline variants...")
    out_filter2 = out_filter1[out_filter1['frequency'] <= 65].copy()
    fig5, ax5, _ = plot_variant_freq_by_sample(out_filter2, sample_order=sample_order)
    fig5.savefig(figdir / "fig5_variant_freq_by_sample_germline_removed.png",  transparent=True, dpi=600)
    plt.close(fig5)

    fig6, ax6, _ = plot_variant_freq_by_sample_gene(
        out_filter2,
        id_column='ID',
        freq_col='frequency',
        gene_col='SYMBOL',
        target_genes=gene_of_interest,
        sample_order=sample_order
    )
    fig6.savefig(figdir / "fig6_variant_freq_by_sample_gene_germline_removed.png",  transparent=True, dpi=600)
    plt.close(fig6)

    # -------------------------------------------------------
    # Apply germline and callability mask to the NGT dataframe 
    # -------------------------------------------------------
    print("Converting variants with a freq >= 65 and callability < 70 in a sample to missing")
    conditions = out_filtered[
        (out_filtered['frequency'] >= 65) | (out_filtered['callability'] < 70)
    ]

    for _, row in conditions.iterrows():
        sample_id = row['sample_id']
        variant_id = row['variant_id']
        if variant_id not in NGT_FILTERED_new.columns:
            continue
        mask = NGT_FILTERED_new['sample_id'].eq(sample_id)
        NGT_FILTERED_new.loc[mask, variant_id] = 3

    vars_invalid = (~NGT_FILTERED_new.iloc[:, :-3].isin([1, 2]).any())
    vars_invalid = (vars_invalid[vars_invalid==True].index.to_list())
    print("Variants entirely removed due to poor callability / germline detection:", len(vars_invalid))

    NGT_FILTERED_new = NGT_FILTERED_new.iloc[:, :-3] 
    if vars_invalid:
        NGT_FILTERED_new.drop(vars_invalid,axis=1, inplace=True)
        DP_new.drop(vars_invalid, axis=1 ,inplace=True)
        AF_FILTERED_new.drop(vars_invalid, axis=1, inplace=True)
        GQ_new.drop(vars_invalid, axis=1, inplace=True)
        annotations.drop(index=vars_invalid, inplace=True)

    #Remove long insertions -> putative false positives (confirmed via IGV visualization)
    length_INS = annotations['ALT'].apply(len)
    long_insertions = length_INS[length_INS > 5].index.to_list()
    print("Long Insertions Removed: ", len(long_insertions))

    if long_insertions:
        NGT_FILTERED_new.drop(long_insertions,axis=1, inplace=True)
        DP_new.drop(long_insertions, axis=1 ,inplace=True)
        AF_FILTERED_new.drop(long_insertions, axis=1, inplace=True)
        GQ_new.drop(long_insertions, axis=1, inplace=True)
        annotations.drop(index=long_insertions, inplace=True)

    # Density plot before binomial thresholding 
    all_af_values = []
    for idx, row in NGT_FILTERED_new.iterrows():
        present_variants = row[row == 1].index
        af_values = AF_FILTERED_new.loc[idx, present_variants].values.flatten()
        all_af_values.extend(af_values)
    all_af_values = np.array(all_af_values)

    fig7, ax = plt.subplots(figsize=(8, 5))
    sns.kdeplot(all_af_values, fill=True, linewidth=2, ax=ax)
    ax.set_xlabel("Allele Frequency (AF)")
    ax.set_ylabel("Density")
    ax.grid(True)
    fig7.savefig(figdir / "fig7_kde_AF_combined_samples_vafref10_vafhet30_vafhom90_postQC.png",  transparent=True, dpi=600, bbox_inches="tight")
    plt.close(fig7)
    print("Density plot before binomial thresholding generated.")

    NGT_final, flips, kept = convert_unbalanced_hets(NGT = NGT_FILTERED_new, DP = DP_new, AF_percent = AF_FILTERED_new, alpha=0.05, p0=0.5, min_dp=10)
    print(f"Flipped: {len(flips)}   Kept: {len(kept)}")
    print("Modified NGT shape:", NGT_final.shape)

    # Denisty plot after binomial thresholding 
    all_af_values = np.array(kept['af_percent'])
    fig8, ax = plt.subplots(figsize=(8, 5))
    sns.kdeplot(all_af_values, fill=True, linewidth=2, ax=ax)
    ax.set_title("Distribution of Allele Frequencies for Heterozygous Variants")
    ax.set_xlabel("Allele Frequency (AF)")
    ax.set_ylabel("Density")
    ax.grid(True)
    fig8.savefig(figdir / "fig8_kde_AF_combined_samples_onesided_binomial_alpha5prct_FDR_vafref10_vafhet30_vafhom90.png",  transparent=True, dpi=600, bbox_inches="tight")
    plt.close(fig8)
    print("Density plot after binomial thresholding generated.")

    vars_binomial_removed = (~NGT_final.isin([1, 2]).any())
    vars_binomial_removed = (vars_binomial_removed[vars_binomial_removed==True].index.to_list())
    print("Variants entirely removed due to binomial thresholding:", len(vars_binomial_removed))

    if vars_binomial_removed:
        NGT_final.drop(vars_binomial_removed,axis=1, inplace=True)
        DP_new.drop(vars_binomial_removed, axis=1 ,inplace=True)
        AF_FILTERED_new.drop(vars_binomial_removed, axis=1, inplace=True)
        GQ_new.drop(vars_binomial_removed, axis=1, inplace=True)
        annotations.drop(index=vars_binomial_removed, inplace=True)

    # Remove variants located within 5 bp of each other in the same cell
    amplicon_id_map = pd.read_csv("/results/rr/study/hg38s/study232-missionbio_TDP-C/outputs/vaf30/amplicon-variant-mapping.tsv", sep="\t", header=0, index_col=0)
    amplicon_id_map = amplicon_id_map.loc[NGT_final.columns]

    NGT_arr= NGT_final.to_numpy() 
    geno_new = mark_close_variants_as_missing(
        geno = NGT_arr,
        variant_ids = NGT_final.columns,
        amplicon_map = amplicon_id_map['amplicon_id'],  
        max_distance = 5
    )
    NGT_final = pd.DataFrame(geno_new, index=NGT_final.index, columns=NGT_final.columns)
    print("Final NGT dimensions:", NGT_final.shape)

    vars_close_remove = (~NGT_final.isin([1, 2]).any())
    vars_close_remove = (vars_close_remove[vars_close_remove==True].index.to_list())
    print("Variants entirely removed due to close distance in the same cell:", len(vars_close_remove))

    if vars_close_remove:
        NGT_final.drop(vars_close_remove,axis=1, inplace=True)
        DP_new.drop(vars_close_remove, axis=1 ,inplace=True)
        AF_FILTERED_new.drop(vars_close_remove, axis=1, inplace=True)
        GQ_new.drop(vars_close_remove, axis=1, inplace=True)
        annotations.drop(index=vars_close_remove, inplace=True)

    # -------------------------------------------------------
    # Final exports
    # -------------------------------------------------------
    print("Exporting final matrices and metadata...")

    # Export files
    print("Checking dimensions of matrices")
    if (NGT_final.shape[0]==AF_FILTERED_new.shape[0]==DP_new.shape[0]==GQ_new.shape[0]==metadata_filtered.shape[0]) & (NGT_final.shape[1]==AF_FILTERED_new.shape[1]==DP_new.shape[1]==GQ_new.shape[1]):
        print("Same set of cells selected across all matrices - Check")
        print("No. of cells = ", NGT_final.shape[0], "No. of variants = ", NGT_final.shape[1])
        print("Exporting matrices...")
        NGT_FILTERED_new.to_parquet(parquet_dir / "NGT_filtered_prebinomial_repeats_germline_freq.parquet",engine="pyarrow",compression="zstd",index=True)
        NGT_final.to_parquet(parquet_dir / "NGT_filtered_final.parquet",engine="pyarrow",compression="zstd",index=True)
        AF_FILTERED_new.to_parquet(parquet_dir / "AF_filtered_final.parquet",engine="pyarrow",compression="zstd",index=True)
        DP_new.to_parquet(parquet_dir / "DP_filtered_final.parquet",engine="pyarrow",compression="zstd",index=True)
        GQ_new.to_parquet(parquet_dir / "GQ_filtered_final.parquet",engine="pyarrow",compression="zstd",index=True)
        metadata_filtered.to_parquet(parquet_dir / "metadata_filtered_repeats_germline_freq.parquet",engine="pyarrow",compression="zstd",index=True)
        flips.to_parquet(parquet_dir / 'flipped_hets_onesided_binomial_alpha5prct_FDR_vafhet30.parquet',engine="pyarrow", compression="zstd", index=True)
        kept.to_parquet(parquet_dir / 'kept_hets_onesided_binomial_alpha5prct_FDR_vafhet30.parquet', engine="pyarrow", compression="zstd", index=True) 
        annotations.to_csv(outdir / 'Annotations_vafhet30.csv', sep='\t')
        print("Modified NGT and reports exported as parquet files.")
    else:
        print("Dimension of matrices do not match. Check for error.")

    # Add back the metadata columns for computing frequencies
    NGT_final['ID'] = metadata_filtered['ID']
    NGT_final['sample_id'] = metadata_filtered['sample_id']
    NGT_final['sample_type'] = metadata_filtered['sample_type']

    ## Plot final frequencies after masking putative germline variants
    out_final = compute_variant_frequencies(NGT_final)
    parts = out_final['variant_id'].astype(str).str.split(':', n=2, expand=True)
    out_final[['chr', 'position', 'var']] = parts
    out_final['position'] = pd.to_numeric(out_final['position'], errors='coerce')

    # Attach SYMBOL to out_final for gene-wise plotting
    gene_lookup = annotations['SYMBOL']
    out_final['SYMBOL'] = out_final['variant_id'].map(gene_lookup)
    out_final = out_final[out_final['frequency'] > 0]

    fig9, ax9, _ = plot_variant_freq_by_sample(out_final, sample_order=sample_order)
    fig9.savefig(figdir / "fig9_variant_freq_targets_by_sample_final.png",  transparent=True, dpi=600)
    plt.close(fig9)

    fig10, ax10, _ = plot_variant_freq_by_sample_gene(
        out_final,
        id_column='ID',
        freq_col='frequency',
        gene_col='SYMBOL',
        target_genes=gene_of_interest,
        sample_order=sample_order
    )
    fig10.savefig(figdir / "fig10_variant_freq_targets_by_gene_final.png",  transparent=True, dpi=600)
    plt.close(fig10)
    
    #Export frequencies
    out_final.to_csv(outdir / "variant_frequencies_filtered_nonzero.tsv",sep="\t")

    #### Extract STG only matrices ####
    NGT_final['sample_tissue'] = metadata_filtered['sample_tissue']
    NGT_final_STG = NGT_final[~NGT_final['sample_tissue'].isin(['FC', 'MTG'])].iloc[:, :-4] 
    
    vars_nonSTG = (~NGT_final_STG.isin([1, 2]).any())
    vars_nonSTG = (vars_nonSTG[vars_nonSTG==True].index.to_list())
    print("Variants entirely removed due to close distance in the same cell:", len(vars_nonSTG))

    NGT_final_STG.drop(vars_nonSTG,axis=1, inplace=True)
    annotations_STG = annotations.drop(index=vars_nonSTG)
    AF_STG = AF_FILTERED_new.loc[NGT_final_STG.index, NGT_final_STG.columns]
    DP_STG = DP_new.loc[NGT_final_STG.index, NGT_final_STG.columns]
    GQ_STG = GQ_new.loc[NGT_final_STG.index, NGT_final_STG.columns]

    if (NGT_final_STG.shape[0]==AF_STG.shape[0]==DP_STG.shape[0]==GQ_STG.shape[0]) & (NGT_final_STG.shape[1]==AF_STG.shape[1]==DP_STG.shape[1]==GQ_STG.shape[1]):
        print("No. of cells = ", NGT_final_STG.shape[0], "No. of variants = ", NGT_final_STG.shape[1])
        print("Exporting matrices...")
    
        NGT_final_STG.to_parquet(parquet_dir / "NGT_filtered_STG_final.parquet",engine="pyarrow",compression="zstd",index=True)
        AF_STG.to_parquet(parquet_dir / "AF_filtered_STG_final.parquet",engine="pyarrow",compression="zstd",index=True)
        DP_STG.to_parquet(parquet_dir / "DP_filtered_STG_final.parquet",engine="pyarrow",compression="zstd",index=True)
        GQ_STG.to_parquet(parquet_dir / "GQ_filtered_STG_final.parquet",engine="pyarrow",compression="zstd",index=True)
        annotations_STG.to_csv(outdir / 'Annotations_STG.csv', sep='\t')

        print("STG matrices exported as parquet files.")
    else:
        print("Dimension of matrices do not match. Check for error.")
    
    log_file.close()

if __name__ == "__main__":
    main()