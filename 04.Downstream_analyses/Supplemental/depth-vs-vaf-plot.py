import missionbio.mosaic as ms
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import matplotlib.gridspec as gridspec

# Plots corresponding to Supplementary Figure S7
def extract_het_values(ngt_matrix, af_matrix, dp_matrix):
    all_vaf = []
    all_depth = []

    for col in range(ngt_matrix.shape[1]):
        ngt_col = ngt_matrix[:, col]
        het_idx = np.where(ngt_col == 1)[0]
        if len(het_idx) >= 1:
            vaf = af_matrix[het_idx, col]
            depth = dp_matrix[het_idx, col]
            all_vaf.extend(vaf)
            all_depth.extend(depth)

    return np.array(all_vaf), np.array(all_depth)


def make_kde_marginal_plot(all_vaf, all_depth, output_svg, output_png):
    fig = plt.figure(figsize=(10, 8))
    gs = gridspec.GridSpec(2, 2, width_ratios=[4, 1], height_ratios=[1, 4], hspace=0.05, wspace=0.05)

    ax_main = fig.add_subplot(gs[1, 0])
    ax_top = fig.add_subplot(gs[0, 0], sharex=ax_main)
    ax_right = fig.add_subplot(gs[1, 1], sharey=ax_main)

    xy = np.vstack([all_vaf, all_depth])
    kde = gaussian_kde(xy)
    x_grid = np.linspace(0, 100, 200)
    y_grid = np.linspace(0, 500, 100)
    X, Y = np.meshgrid(x_grid, y_grid)
    Z = kde(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)

    ax_main.contourf(X, Y, Z, levels=20, cmap='Reds')
    ax_main.contour(X, Y, Z, levels=20, colors='darkred', linewidths=0.3, alpha=0.5)
    ax_main.axvline(x=30, color='red', linestyle='dotted', linewidth=2)
    ax_main.axhline(y=10, color='red', linestyle='dotted', linewidth=2)
    ax_main.set_xlim(-2, 100)
    ax_main.set_ylim(-30, 500)
    ax_main.set_xlabel('Allelic Fraction', fontsize=20)
    ax_main.set_ylabel('Depth', fontsize=20)
    ax_main.set_xticks(range(0, 101, 20))
    ax_main.set_yticks(range(0, 501, 100))
    ax_main.spines['right'].set_visible(False)
    ax_main.spines['top'].set_visible(False)
    for spine in ax_main.spines.values():
        spine.set_linewidth(2)
    ax_main.tick_params(axis='both', labelsize=18, width=2)

    kde_vaf = gaussian_kde(all_vaf)
    x_range = np.linspace(0, 100, 300)
    ax_top.fill_between(x_range, kde_vaf(x_range), alpha=0.4, color='darkred')
    ax_top.plot(x_range, kde_vaf(x_range), color='darkred', linewidth=1.5)
    ax_top.axvline(x=30, color='red', linestyle='dotted', linewidth=2)
    ax_top.set_xlim(-2, 100)
    ax_top.spines['right'].set_visible(False)
    ax_top.spines['top'].set_visible(False)
    ax_top.spines['left'].set_visible(False)
    ax_top.tick_params(labelbottom=False, left=False, labelleft=False, width=2)

    kde_depth = gaussian_kde(all_depth)
    y_range = np.linspace(0, 500, 300)
    ax_right.fill_betweenx(y_range, kde_depth(y_range), alpha=0.4, color='darkred')
    ax_right.plot(kde_depth(y_range), y_range, color='darkred', linewidth=1.5)
    ax_right.axhline(y=10, color='red', linestyle='dotted', linewidth=2)
    ax_right.set_ylim(-30, 500)
    ax_right.spines['right'].set_visible(False)
    ax_right.spines['top'].set_visible(False)
    ax_right.spines['bottom'].set_visible(False)
    ax_right.tick_params(labelleft=False, bottom=False, labelbottom=False, width=2)

    fig.savefig(output_svg, bbox_inches='tight', transparent=True)
    fig.savefig(output_png, dpi=600, bbox_inches='tight', transparent=True)
    plt.close(fig)


# Load pre-QC data from the integrated h5 file
file = "/home/AD/vbidhan/study232-missionbio_TDP-C/data/outputs/integrated_h5_files/hg38/combined_samples_allvars_vafref10_vafhet30_vafhom90.h5"
data_preQC = ms.load(file, raw=False, filter_variants=False, filter_cells=False, single=True)

# Load filtered data from parquet files (post QC)
parquet_dir = "/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/parquets/"
AF_final_df = pd.read_parquet(f"{parquet_dir}AF_filtered_final.parquet")
DP_final_df = pd.read_parquet(f"{parquet_dir}DP_filtered_final.parquet")
NGT_final_df = pd.read_parquet(f"{parquet_dir}NGT_filtered_final.parquet")

print(NGT_final_df.shape)
print(DP_final_df.shape)
print(AF_final_df.shape)

# Prepare the pre-QC matrices (removing the SNCA region)
NGT_prefilter = data_preQC.dna.get_attribute('NGT')
col_names = list(NGT_prefilter)
snca_idx_remove = [
    i for i, col in enumerate(col_names)
    if col.startswith("chr4:") and 89726454 <= int(col.split(":")[1]) <= 89835667
]
rows_to_remove = NGT_prefilter.index[~NGT_prefilter.index.isin(NGT_final_df.index)].tolist()
row_idx = [NGT_prefilter.index.get_loc(i) for i in rows_to_remove]

NGT_preQC = data_preQC.dna.layers['NGT']
AF_preQC = data_preQC.dna.layers['AF']
DP_preQC = data_preQC.dna.layers['DP']

NGT_preQC = np.delete(NGT_preQC, snca_idx_remove, axis=1)
NGT_preQC = np.delete(NGT_preQC, row_idx, axis=0)
AF_preQC = np.delete(AF_preQC, snca_idx_remove, axis=1)
AF_preQC = np.delete(AF_preQC, row_idx, axis=0)
DP_preQC = np.delete(DP_preQC, snca_idx_remove, axis=1)
DP_preQC = np.delete(DP_preQC, row_idx, axis=0)

print(NGT_preQC.shape)
print(AF_preQC.shape)
print(DP_preQC.shape)

# Convert filtered dataframes to numpy arrays
NGT_final = NGT_final_df.to_numpy()
DP_final = DP_final_df.to_numpy()
AF_final = AF_final_df.to_numpy()

# Extract heterozygous values for pre-QC data and make the corresponding plot
all_vaf_preQC, all_depth_preQC = extract_het_values(NGT_preQC, AF_preQC, DP_preQC)
print('All VAF dimensions pre QC:', all_vaf_preQC.shape)
make_kde_marginal_plot(
    all_vaf_preQC,
    all_depth_preQC,
    "/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/figures/supplementary/vaf_plots/2dkde_marginal_vaf_depth.svg",
    "/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/figures/supplementary/vaf_plots/2dkde_marginal_vaf_depth.png"
)

# Extract heterozygous values for post-QC data and make the corresponding plot
all_vaf_postQC, all_depth_postQC = extract_het_values(NGT_final, AF_final, DP_final)
print('All VAF dimensions post QC:', all_vaf_postQC.shape)
make_kde_marginal_plot(
    all_vaf_postQC,
    all_depth_postQC,
    "/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/figures/supplementary/vaf_plots/2dkde_marginal_vaf_depth_postQC.svg",
    "/home/AD/vbidhan/study232-missionbio_TDP-C/manuscript_data/figures/supplementary/vaf_plots/2dkde_marginal_vaf_depth_postQC.png"
)


