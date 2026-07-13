# Import mosaic libraries
import missionbio.mosaic as ms
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.spatial.distance import cdist, pdist, euclidean
from matplotlib.patches import Circle
import itertools
from pathlib import Path
import os
import sys
import colorsys
from collections import defaultdict


def convert_to_int(value):
    if value.name == 'ID':
        return value
    return value.astype(np.int8)

def count_threes(axis):
    count = 0
    for value in axis:
        if value == 3:
            count += 1
    return count

# Load demultiplexing variants
def extract_demux_variants(vcf_path):
    variants = []
    vcf_path = Path(vcf_path)

    with vcf_path.open() as f:
        for line in f:
            if line.startswith("#"):
                continue

            chrom, pos, _id, ref, alt, *_ = line.rstrip("\n").split("\t")
            for a in alt.split(","):
                variants.append(f"{chrom}:{pos}:{ref}/{a}")

    return variants

def replace_with_conditions(column):
    mean_value = column.mean()
    if mean_value > 0 and mean_value < 0.5:
        return 0
    elif mean_value >= 0.5 and mean_value < 1.5:
        return 1
    elif mean_value >= 1.5:
        return 2
    else:
        return mean_value  

def select_n_components(explained_var, min_pcs=5, target=0.80):
    cumvar = np.cumsum(explained_var)
    n_target = np.searchsorted(cumvar, target) + 1
    return max(min_pcs, n_target)
    
def sample_two_cells(df, num_samples=2):
    while True:
        sample = df.sample(n=num_samples, replace=False)
        first_chars = [col[0] for col in sample.T.columns]
        if len(set(first_chars)) == len(first_chars):  
            return sample

def get_n_clusters(poolNumber):
    return {
        10: 8,
        8: 7,
        5: 7,
        6: 5,
        11: 5,
        14: 5, 
        15: 5, 
        20: 5
    }.get(poolNumber, 6)


# Get minimum number of representative cells to be sampled 
def get_min_ncell(p_cells, min_size=5, max_cap=50):
    if not p_cells:
        raise ValueError("p_cells is empty")

    lengths = sorted(len(c) for c in p_cells)

    for l in lengths:
        if l >= min_size:
            return min(l, max_cap)

    return min(lengths[-1], max_cap)

class Tee:
    def __init__(self, *files):
        self.files = files

    def write(self, data):
        for f in self.files:
            f.write(data)
            f.flush()

    def flush(self):
        for f in self.files:
            f.flush()


demux_vars = extract_demux_variants("/home/vanshika/study232-missionbio_TDP-C/data/final/hg38/demux_vars_hg38.vcf")
output_dir_path = "/home/vanshika/study232-missionbio_TDP-C/data/final/hg38/demultiplexing_results/"
pools = [4,5,6,7,8,9,10,11,12,13,14,15,16,17,20,22,23,24,25, "Test"]
SEED = 10
np.random.seed(SEED)

for pool_number in pools:
    print(f"Analysing pool{pool_number}")
    output_dir= Path(os.path.join(output_dir_path, f"Pool{pool_number}/"))
    output_dir.mkdir(parents=True, exist_ok=True)

    original_stdout = sys.stdout
    original_stderr = sys.stderr

    log_file = Path(output_dir) / "qc_output_dnavarsonly.log"
    with open(log_file, "w") as log:
        sys.stdout = Tee(original_stdout, log)
        sys.stderr = Tee(original_stderr, log)

        output_barcodes = Path(os.path.join(output_dir_path, "barcodes"))
        output_barcodes.mkdir(parents=True, exist_ok=True)

        ###Set theresholds
        #SNPs which are genotyped for x% cells
        if (pool_number==24):
            demux_coverage = 90 
        else:
            demux_coverage = 85 

        n_clusters = get_n_clusters(pool_number)
        print(f"Number of clusters: {n_clusters}")
        THRESH = 0.10 # Threshold for % cells mapping to a cluster 
        
        # Minimum number of confident cells required to assign a sample to a cluster
        if pool_number == 6:
            min_confident_cells = 6 #To exclude sample C2_C
        else:
            min_confident_cells = 5
        
        h5path = f'/home/AD/vbidhan/study343-tmp_vbidhan/tapestri_v2_runs/hg38/pool{pool_number}/results/pool{pool_number}-v2.dna.h5'
        sample = ms.load(h5path, raw=False, filter_variants=False, filter_cells=False) 
        cells_preQC = sample.dna.shape[0]
        print("Number of cells:", sample.dna.shape[0], ", Number of variants:", sample.dna.shape[1])
        
        # Select appropriate demultiplexing SNPs for clustering
        dds = sample.dna[:, demux_vars]
        dna_vars = dds.filter_variants(min_dp=10,min_gq=30,min_prct_cells=80)
        print(f"Number of demultiplexing variants after filtering: {len(dna_vars)}")
        dds.filter_variants( min_dp=10,min_gq=30,min_prct_cells=80)
        dds = dds[:, dna_vars]
        
        # Find unreliable cells with high missingness as they cannot be classified reliably and remove them from the analysis
        data = dds.layers['NGT_FILTERED']
        NGT = pd.DataFrame(data)
        NGT.columns = dds.col_attrs['id']
        NGT.index = dds.row_attrs['barcode']
        cells_missing_data = pd.DataFrame((NGT == 3).sum(axis=1))
        cells_missing_data.columns = ["missingness"] 
        unreliable_cells = cells_missing_data[cells_missing_data['missingness'] > 10].index #excludes cells with more than 10 missing variants 
        
        if (len(unreliable_cells) > 0):
            sample.dna = sample.dna.drop(unreliable_cells)
            NGT.drop(unreliable_cells, axis=0, inplace=True)
        
        #Export the list of unreliable cells (not used for clustering downstream)
        with open(f"{output_dir}/unreliable_cells.txt", "w") as f:
            for item in unreliable_cells:
                f.write(f"{item}\n")

        path = f"/home/AD/vbidhan/study232-missionbio_TDP-C/data/genotypes-per-pool/hg38/pool{pool_number}.csv"
        variants = pd.read_csv(path)
        variants.set_index('ID', inplace=True)
        if (pool_number == 9):
            variants.drop(['chr10:111434641:C/G', 'chr17:43721456:G/T'], axis = 0, inplace=True) #These SNPs do not match the consensus genotypes for cluster F3_P
        
        # Samples without germline SNPs known
        fixes = {
            7:  ('D1_C', 5)
        }
        if pool_number in fixes:
            sample_name, _ = fixes[pool_number]
            variants.drop(columns=[sample_name], inplace=True, errors="ignore")

        variants = variants.loc[[v for v in demux_vars if v in variants.index]]
        missingness_per_SNP = variants.apply(count_threes, axis=1)
        
        # Select SNPs for which genotypes are available for all samples
        covered_snps = missingness_per_SNP[missingness_per_SNP == 0].index.to_list()
        
        # Find true sample representative cells
        covered_snps_object = sample.dna[: , covered_snps]
        covered_snps_object = covered_snps_object[:, covered_snps_object.filter_variants(min_prct_cells=demux_coverage)] 
        NGT_covered_snps = pd.DataFrame(covered_snps_object.layers['NGT'])
        NGT_covered_snps.columns = covered_snps_object.col_attrs['id']
        NGT_covered_snps.index = covered_snps_object.row_attrs['barcode']

        demux_data = variants.T[covered_snps]
        demux_data = demux_data[NGT_covered_snps.columns] # subset for SNPs genotyped for most of the cells 
        demux_data = demux_data.reindex(columns=NGT_covered_snps.columns)
        demultiplexing_SNPs = demux_data.columns
        demux_data['GT_code'] = demux_data.apply(lambda row: ''.join(map(str, row)), axis=1)
        NGT_covered_snps['code'] = NGT_covered_snps.apply(lambda row: ''.join(map(str, row)), axis=1)
        print("Number of variants used for demultiplexing:" f"{len(demultiplexing_SNPs)}")

        known_gts_cells = []
        print("\nNumber of sample representative cells identified:")
        for i in range(len(demux_data.index)):
            variant_code = str(demux_data.iloc[i, -1])
            perfect_cells = NGT_covered_snps[NGT_covered_snps['code'] == variant_code].index
            known_gts_cells.append(perfect_cells)
            print(f"{demux_data.index[i]} = ",  len(perfect_cells))
        print("\n")

        p_cells = list(known_gts_cells)

        #Mean impute the missing genotypes for demux variants prior to PCA and clustering
        NGT.replace(3, np.nan, inplace=True)
        modified_means = NGT.apply(replace_with_conditions)
        NGT.fillna(modified_means, inplace=True)

        pca = PCA(n_components=len(dna_vars))
        pca_data = pca.fit_transform(NGT)
        pca_data = pd.DataFrame(pca_data)
        pca_data.index = sample.dna.row_attrs['barcode']

        explained_variance = pca.explained_variance_ratio_

        n_comp = select_n_components(explained_variance)
        print("No. of components selected:", n_comp)

        pca = PCA(n_components=n_comp, svd_solver="randomized", random_state=SEED)
        pca_data = pca.fit_transform(NGT)
        pca_data = pd.DataFrame(pca_data)
        pca_data.index = sample.dna.row_attrs['barcode']

        ## PCA Plot 
        plt.figure(figsize=(8, 6))
        scatter = plt.scatter(pca_data[0], pca_data[1], edgecolors='k', s=50, c='#8C8C8C')
        plt.title(f'PCA - Based on {n_comp} components')
        plt.xlabel('PCA axis 1')
        plt.ylabel('PCA axis 2')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/1_PCA_initial.png", transparent=True, dpi=600)
        plt.close()

        colors = plt.cm.tab10.colors[:(n_clusters + 2)]
        ## PCA Plot to highlight sample representative cells
        plt.figure(figsize=(8, 6))
        plt.scatter(pca_data[0], pca_data[1], edgecolors='k',s=50,c='#8C8C8C')
        legend_handles = []
        for i, cells in enumerate(p_cells):
            color = colors[i]
            # Highlight points
            plt.scatter(pca_data.loc[cells, 0],pca_data.loc[cells, 1],edgecolors='k',c=[color],s=50,marker='o')
            legend_handles.append(Circle((0.5, 0.5), 0.2, color=color, label=demux_data.index[i]))

        legend = plt.legend(handles=legend_handles, title='Samples',  loc='upper right', frameon=False)
        legend.get_frame().set_alpha(0.0)
        legend.get_frame().set_edgecolor('none')
        plt.title(f'Pool {pool_number}- Sample representative cells')
        plt.xlabel('PCA axis 1')
        plt.ylabel('PCA axis 2')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/2_PCA_sample_representative_cells.png", transparent=True, dpi=600)
        plt.close()

        # Running k-means algorithm 
        kmeans = KMeans(n_clusters=n_clusters, random_state=10, n_init = 30) 
        kmeans.fit(NGT)
        labels = kmeans.labels_
        cluster_assignments= np.array(labels, dtype=float).astype(int)
        cluster_assignments = cluster_assignments.tolist()
        #Add labels to the main "dds" dataframe
        dds.row_attrs['label'] = labels

        cluster_assignments_df = pd.DataFrame(cluster_assignments)
        cluster_assignments_df.index = NGT.index

        cluster_cells = {
            i: cluster_assignments_df[cluster_assignments_df[0] == i].index.tolist()
            for i in range(n_clusters)
        }
    
        NGT['Cluster'] = dds.row_attrs['label']
        
        print(f"\nSample -> Cluster mappings:")
        # Sample to Cluster mappings
        sample_to_cluster  = {}   # sample -> main cluster to which it is mapped
        sample_to_clusters = {}   # sample -> list of clusters to which it is mapped
        sample_to_props    = {}   # sample -> {cluster: proportion} (full, unfiltered)

        for i, cells in enumerate(p_cells):
            s = demux_data.index[i]
            if (len(cells) >= min_confident_cells):
                cluster_props = NGT.loc[cells, "Cluster"].value_counts(normalize=True)
                sample_to_props[s] = cluster_props.to_dict()

                if cluster_props.empty:
                    print(f"Sample {s}: no confident cells found")
                    sample_to_cluster[s] = None
                    sample_to_clusters[s] = []
                    continue

                # majority cluster always kept
                top_cluster = cluster_props.idxmax()
                sample_to_cluster[s] = top_cluster

                # keep additional clusters only if they are >= 10%
                kept = cluster_props[cluster_props >= THRESH].index.tolist()

                # ensure majority is always included even if < THRESH (edge case)
                if top_cluster not in kept:
                    kept = [top_cluster] + kept

                # optional: sort kept clusters by proportion descending
                kept = sorted(kept, key=lambda c: cluster_props.loc[c], reverse=True)

                sample_to_clusters[s] = kept

                if len(kept) == 1:
                    print(f"Sample {s} maps to cluster {top_cluster}")
                else:
                    print(f"Sample {s} maps to multiple clusters (>= {THRESH*100:.0f}%):")
                    for c in kept:
                        print(f"  Cluster {c}: {cluster_props.loc[c]*100:.2f}%")
            else:
                print(f"Sample {s} has less than {min_confident_cells} confident cells; skipping assignment")        
        

        print(f"\nCluster -> Sample mappings:")
        # Cluster to Sample mappings
        all_clusters = sorted(NGT["Cluster"].dropna().unique())
        cluster_to_samples  = {c: {} for c in all_clusters}   # cluster -> {sample: proportion}
        cluster_to_majority = {c: [] for c in all_clusters}   # cluster -> [samples where this is majority]

        for i, cells in enumerate(p_cells):
            s = demux_data.index[i]
            if (len(cells) >= min_confident_cells):
                cluster_props = NGT.loc[cells, "Cluster"].value_counts(normalize=True)

                if cluster_props.empty:
                    continue

                top_cluster = cluster_props.idxmax()
                cluster_to_majority[top_cluster].append(s)

                # only include sample under a cluster if that cluster constitutes ≥ 10% cells of the sample
                for cluster, prop in cluster_props.items():
                    if prop >= THRESH:
                        cluster_to_samples[cluster][s] = prop
            else:
                print(f"Note: Sample {s} has less than {min_confident_cells} confident cells; skipping assignment")

        for cluster in all_clusters:
            if not cluster_to_samples[cluster]:
                print(f"Cluster {cluster}: No samples mapped")
                continue
            for s, prop in sorted(
                cluster_to_samples[cluster].items(),
                key=lambda x: x[1],
                reverse=True
            ):
                print(f"Cluster {cluster}: Sample {s}- {prop*100:.2f}%")
        print("\n")

        samples = list(demux_data.index) 
        if pool_number in fixes:
            s, cluster = fixes[pool_number]
            sample_to_cluster[s] = cluster
            sample_to_clusters[s] = [cluster]
            cluster_to_majority[cluster] = [s]
            cluster_to_samples[cluster] = {s: 1.0}
            samples.append(s)

        # Select colors for plotting 
        doublet_color = colors[-1]
        sample_colors = colors[:-1] 
        cluster_colors = [doublet_color] * n_clusters  # default every cluster to doublet cluster's color 
        for i, s in enumerate(samples):
            clusters = sample_to_clusters.get(s, [])  # <-- list of clusters (multi-mapped)
            if not clusters:
                continue
            if i >= len(sample_colors):
                raise ValueError("Not enough sample_colors (need <= len(samples) real samples).")
            for cl in clusters:
                cl = int(cl)
                if cluster_colors[cl] != doublet_color:
                    continue
                cluster_colors[cl] = sample_colors[i]
        
        # Color map with selected colors
        cmap = ListedColormap(cluster_colors)

        # PCA coloured according to k-means clustering
        plt.figure(figsize=(8, 6))
        scatter = plt.scatter(pca_data[0], pca_data[1], edgecolors='k', s=50, c=cluster_assignments, cmap=cmap)
        legend_labels = sorted(set(cluster_assignments))
        legend = plt.legend(handles=scatter.legend_elements()[0], labels=legend_labels, title='Clusters',  loc='upper right', frameon=False)

        legend.get_frame().set_alpha(0.0)
        legend.get_frame().set_edgecolor('none')

        plt.title(f'Pool {pool_number} - Kmeans Clustering')
        plt.xlabel('PCA axis 1')
        plt.ylabel('PCA axis 2')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/3_PCA_kmeans_clustering.png", transparent=True, dpi=600)
        plt.close()

        #Doublet Detection
        cluster_ncell = NGT.groupby('Cluster').size().reset_index(name='n_cells')

        # All clusters that appear in any sample mapping (multi-mapped)
        all_mapped_clusters = {
            int(cl)
            for clusters in sample_to_clusters.values()
            for cl in clusters
        }

        # Clusters that are majority in any sample
        mapped_majority_clusters = {v for v in sample_to_cluster.values() if v is not None}

        # Clusters that do not map to any sample
        missing_clusters = sorted(set(range(n_clusters)) - all_mapped_clusters)

        # A list of cells that do not match any sample (doublets or otherwise)
        doublets = []
        for i in range(len(missing_clusters)):
            cluster = missing_clusters[i]
            doublets.extend(cluster_cells[cluster])

        # Select confident cells from the samples that are matched to clusters
        selected_known_gts_cells = [
            cells
            for label, cells in zip(demux_data.index, known_gts_cells)
            if label in sample_to_cluster
        ]
        p_cells_filt = list(selected_known_gts_cells)
        min_ncell = get_min_ncell(p_cells_filt) 
        print("Minimum number of representative cells sampled:", min_ncell)
        
        # List of all representative cells
        all_representative_cells = [cell for cells in p_cells_filt for cell in cells]

        if pool_number in fixes:
            _, cluster = fixes[pool_number]
            mapped_majority_clusters = mapped_majority_clusters - {cluster}
        
        sampled_pool_list = {
            cluster_name: group[group.index.isin(all_representative_cells)].sample(min_ncell, replace=False)
            for cluster_name, group in NGT.groupby("Cluster")
            if cluster_name in mapped_majority_clusters #excudes minor clusters that might be mapped to the same sample (so a sample is only represented once - avoids creating doublets between two same samples)
        }

        modified_dfs = []
        for idx, df in sampled_pool_list.items():
            df.index = [f"{idx}_{row_name}" for row_name in df.index.astype(str)]
            modified_dfs.append(df)

        sampled_pool = pd.concat(modified_dfs)
        sampled_pool.drop('Cluster', axis =1, inplace =True)

        ## Generate artificial doublets (10%)
        total_doublets = int(len(NGT) * THRESH)
        all_doublets = []

        for i in range(total_doublets):
            sampled_cells = sample_two_cells(sampled_pool)
            mean_col = sampled_cells.mean(axis=0) 
            df= pd.DataFrame(mean_col)
            new_col_name = '_'.join(sampled_cells.T.columns)
            df.columns = [new_col_name]
            all_doublets.append(df)
            
        artificial_doublets_df = pd.concat(all_doublets, axis = 1, ignore_index=False).T

        transformed_new_points = pca.transform(artificial_doublets_df)
        transformed_new_points = pd.DataFrame(transformed_new_points)
        transformed_new_points.index = artificial_doublets_df.index

        combinations = list(itertools.combinations(list(mapped_majority_clusters), 2))
        no_artificial_clusters = len(combinations)
        print("Number of artificial doublet clusters:", no_artificial_clusters)

        #Perform k-means clustering on artificial doublets
        kmeans_sampled = KMeans(n_clusters=no_artificial_clusters, n_init=30, random_state=10)
        kmeans_sampled.fit(artificial_doublets_df)
        cluster_assignments_artificial_cells = kmeans_sampled.labels_
        cluster_centers_sampled = kmeans_sampled.cluster_centers_
        pca_data['Cluster'] = dds.row_attrs['label']

        #Visualise artificial doublets in PCA space
        artifical_cluster_colours = [
            colorsys.hsv_to_rgb(h, 0.7, 0.9)
            for h in np.linspace(0, 1, 12)][-no_artificial_clusters:]  # take 10, avoiding tab10-like hues

        plt.figure(figsize=(8, 6))
        cmap_artificial = ListedColormap(artifical_cluster_colours)
        plt.scatter(pca_data[0], pca_data[1], edgecolors='k', s=50, color='grey')
        scatter_artificial_doublets = plt.scatter(transformed_new_points[0], transformed_new_points[1], label='New Points', c=cluster_assignments_artificial_cells, edgecolors='k', cmap=cmap_artificial)
        legend_labels = sorted(set(cluster_assignments_artificial_cells))
        labels = [str(i) for i in legend_labels]
        legend = plt.legend(handles=scatter_artificial_doublets.legend_elements()[0], labels=labels, title='Clusters', loc='upper right')
        legend.get_frame().set_alpha(0.0)
        legend.get_frame().set_edgecolor('none')
        plt.title(f'Pool {pool_number} - Artificial Cell clusters')
        plt.xlabel('PCA axis 1')
        plt.ylabel('PCA axis 2')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/4_PCA_artificial_clusters.png", transparent=True, dpi=600)
        plt.close()

        # Centeroids of pca data
        cluster_centers = pca.transform(kmeans.cluster_centers_)
        Cluster_centers = []
        for i in range(cluster_centers.shape[0]):
            center = cluster_centers[i].reshape(1, -1)
            Cluster_centers.append(center)

        # Calculate distances from real cells to their cluster centroids
        distances_real = []
        for cluster in range(n_clusters):
            dist_real = pd.DataFrame(cdist(pca_data[pca_data['Cluster']==cluster].iloc[:,:-1], Cluster_centers[cluster], metric='euclidean'))
            dist_real.index = pd.DataFrame(pca_data[pca_data['Cluster']==cluster]).index  
            #dist_real['Cluster'] = cluster
            dist_real.columns = ["distance"]
            distances_real.append(dist_real)
        
        # Calculate distances from real cells to artificial cluster centroids
        distances_artificial = {i: {j: [] for j in range(no_artificial_clusters)} for i in range(n_clusters)}
        cluster_centers_reduced = pca.transform(cluster_centers_sampled)
        for i in range(n_clusters):
            for j in range(no_artificial_clusters):
                #Check distance of cluster 0 cells from artificial cluster centers
                dist = pd.DataFrame(cdist(pca_data[pca_data['Cluster']==i].iloc[:,:-1], [cluster_centers_reduced[j]], metric='euclidean'))
                dist.index = pd.DataFrame(pca_data[pca_data['Cluster']==i]).index 
                dist.columns = ["distance"]
                dist['real_cluster'] = i
                dist['artificial_cluster'] = j
                distances_artificial[i][j].append(dist)
        
        # Make a dataframe with distance of each cell to the center of the artificial cell cluster
        cluster_dfs = {}
        for cluster in list(mapped_majority_clusters):  #cluster already assumed to be doublets is not used 
            dfs = []
            for i in range(no_artificial_clusters):
                df = pd.DataFrame(distances_artificial[cluster][i][0]['distance'])
                df.columns = [f"distance_{cluster}{i}"]
                dfs.append(df)
            distances_cluster_artificial = dfs[0]
            for df in dfs[1:]:
                distances_cluster_artificial = distances_cluster_artificial.join(df, how='outer')
            cluster_dfs[f'cluster{cluster}_artificial'] = distances_cluster_artificial

        # Identify outliers based on distances to artificial cluster centroids
        outliers_highCF = []
        outliers_medCF = []
        outliers_lowCF = []

        for cluster_index in list(mapped_majority_clusters):
            cluster_key = f"cluster{cluster_index}_artificial"
            
            data = cluster_dfs[cluster_key]
            
            for col_index in range(len(data.columns)):
                outliers_col_1 = data[data.iloc[:, col_index] < distances_real[cluster_index]['distance']].index.to_list()
                outliers_col_2 = data[data.iloc[:, col_index] < 2*distances_real[cluster_index]['distance']].index.to_list()
                outliers_col_3 = data[data.iloc[:, col_index] < 3*distances_real[cluster_index]['distance']].index.to_list()
                outliers_highCF.extend(outliers_col_1)
                outliers_medCF.extend(outliers_col_2)
                outliers_lowCF.extend(outliers_col_3)

        outliers_highCF = list(set(outliers_highCF))
        outliers_medCF = list(set(outliers_medCF) - set(outliers_highCF))
        outliers_lowCF = list(set(outliers_lowCF) - set(outliers_medCF) - set(outliers_highCF))

        print(f"Number of high CF outliers: {len(outliers_highCF)}")
        print(f"Number of doubelts/unmapped cells: {len(doublets)}")

        # Visualise doublet cells in PCA space
        plt.figure(figsize=(8, 6))
        plt.scatter(pca_data[0], pca_data[1], edgecolors='k', s=50, color='grey')
        plt.scatter(pca_data.loc[outliers_highCF , 0], pca_data.loc[outliers_highCF , 1], edgecolors='k', s=50, color='red')
        plt.scatter(pca_data.loc[doublets , 0], pca_data.loc[doublets, 1], edgecolors='k', s=50, color='red')
        plt.title(f'Doublets- Pool {pool_number}')
        plt.xlabel('PCA axis 1')
        plt.ylabel('PCA axis 2')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/5_PCA_highCF_doublets.png", transparent=True, dpi=600)
        plt.close()

        # Add cell identity as metadata
        metadata = pd.DataFrame(NGT['Cluster'])
        metadata['sample_id'] = "unassigned"
        metadata['doublet_pred'] = "Singlet"

        cluster_to_samples = defaultdict(list)
        for s, clusters in sample_to_clusters.items():
            for cl in clusters:
                cluster_to_samples[int(cl)].append(s)

        # Cluster -> label. If multiple samples map to a cluster, join them.
        cluster_to_label = {
            cl: "|".join(sorted(set(samples)))   # "SampleA|SampleB"
            for cl, samples in cluster_to_samples.items()
        }

        # assign sample_id by cluster
        metadata["sample_id"] = metadata["Cluster"].map(cluster_to_label)
        # mark clusters that map to nothing as Doublet
        metadata.loc[metadata["Cluster"].isin(missing_clusters), "sample_id"] = "Doublet"

        # Define entries in 'doublet_pred' column.
        metadata.loc[metadata['Cluster'].isin(missing_clusters),'doublet_pred'] = "High"
        metadata.loc[outliers_highCF, 'doublet_pred'] = "High"
        metadata.loc[outliers_medCF, 'doublet_pred'] = "Medium"
        metadata.loc[outliers_lowCF, 'doublet_pred'] = "Low"

        metadata.to_csv(f"{output_barcodes}/pool{pool_number}.csv", index = True)

        cells_to_remove = metadata[metadata['doublet_pred'].isin(['High'])].index.tolist() 
        doublet_rate = (len(cells_to_remove)/metadata.shape[0] * 100)
        metadata_cleaned = metadata.drop(cells_to_remove, axis = 0)
        metrics = metadata_cleaned.groupby('sample_id').size()

        # Sample- level metrics
        rows = []
        for i in range(len(demux_data)):
            sample_id = demux_data.index[i]
            if sample_id in metrics.index:
                total_cells = metrics.loc[sample_id]
            else:
                total_cells = 0
            removed_count = sum(cell in cells_to_remove for cell in p_cells[i])
            true_reps = len(known_gts_cells[i])
            rows.append({
                "Pool": pool_number,
                "Sample_id": sample_id,
                "Total_cells": total_cells,
                "True_cells": true_reps,
                "False_doublets": removed_count,
                "False_doublets_rate_prct": (
                    np.round(removed_count / true_reps * 100, 2)
                    if true_reps > 0 else 0.0
                )
            })
        sample_demux_info = pd.DataFrame(rows).set_index("Pool")

        if pool_number in fixes:
            sample_name, _ = fixes[pool_number]
            entry = {
            'Sample_id': sample_name,
            'Total_cells': metrics.loc[sample_name],
            }
            sample_demux_info = pd.concat([sample_demux_info, pd.DataFrame([entry], index = [pool_number])], axis=0)    
        sample_demux_info.to_csv(f"{output_dir}/sample_demux_metrics.csv", index=True)

        # Print pool metrics
        print("\nPool-level metrics:")
        total_cells_removed = len(cells_to_remove) + len(unreliable_cells)
        print("No. of cells pre QC:", cells_preQC)
        print("No. of cells post QC:", metadata_cleaned.shape[0])
        print("No. of unreliable cells removed:", len(unreliable_cells))
        print("No. of cells removed as doublets:", len(cells_to_remove))
        print("Total no. of cells removed:", total_cells_removed)
        print("Doublet rate (%):", np.round(doublet_rate, 2))
        print("Cells (%) lost during QC:", np.round((total_cells_removed/cells_preQC * 100), 2))

        for i in range(len(p_cells)):
            p_cells[i] = [cell for cell in p_cells[i] if cell not in cells_to_remove]
        NGT.drop(cells_to_remove, axis=0, inplace=True)
        cluster_assignments_df.drop(cells_to_remove, axis=0, inplace=True)
        pca_data.drop(cells_to_remove, axis=0, inplace=True)
        sample.dna = sample.dna.drop(cells_to_remove)

        # Colour mapping for the cleaned data
        NGT["sample_id"] = NGT.index.map(metadata_cleaned["sample_id"])

        label_order = NGT['sample_id'].unique().tolist()
        ordered_colors = []
        for label in label_order:
            idx = sample_to_cluster.get(label)
            if idx is None:
                continue  # or append a default color if you prefer
            ordered_colors.append(cluster_colors[idx])

        labels_new = NGT["sample_id"]
        cluster_assignments_new = np.asarray(labels_new)
        codes = pd.Categorical(labels_new,categories=label_order,ordered=True).codes
        cmap_ordered = ListedColormap(ordered_colors)

        ## PCA after QC
        plt.figure(figsize=(8, 6))
        scatter = plt.scatter(pca_data[0], pca_data[1], edgecolors='k', s=50, c=codes, cmap=cmap_ordered)
        legend = plt.legend(handles = scatter.legend_elements()[0], labels=label_order, title="Sample ID", loc='upper right', frameon=False)
        legend.get_frame().set_alpha(0.0)
        legend.get_frame().set_edgecolor('none')
        plt.title(f'Pool {pool_number}- Post QC')
        plt.xlabel('PCA axis 1')
        plt.ylabel('PCA axis 2')        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/6_PCA_final.png", transparent=True, dpi=600)
        plt.close()

        ## UMAP after QC
        # UMAP
        dds.find_clones(similarity=1)
        umap_data = pd.DataFrame(dds.row_attrs['umap'])
        umap_data.index = dds.row_attrs['barcode']
        umap_data = umap_data.loc[pca_data.index] 

        plt.figure(figsize=(8, 6))
        scatter = plt.scatter(umap_data[0], umap_data[1], edgecolors='k', s=50, c=codes, cmap=cmap_ordered)
        legend = plt.legend(handles = scatter.legend_elements()[0], labels=label_order, title="Sample ID", loc='upper right', frameon=False)
        legend.get_frame().set_alpha(0.0)
        legend.get_frame().set_edgecolor('none')
        plt.title(f'Pool {pool_number}- Post QC UMAP')
        plt.xlabel('UMAP axis 1')
        plt.ylabel('UMAP axis 2')        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/7_UMAP_final.png", transparent=True, dpi=600)
        plt.close()

    sys.stdout = original_stdout
    sys.stderr = original_stderr

    print(f"Done pool {pool_number}. Log: {log_file}")
    print("\n")
