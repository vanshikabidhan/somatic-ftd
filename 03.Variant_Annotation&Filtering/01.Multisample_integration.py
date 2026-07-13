import missionbio.mosaic as ms
import numpy as np
import pandas as pd
import sys
print("Packages loaded")

def setup_logging(log_path: str):
    log_file = open(log_path, 'w')
    sys.stdout = log_file
    return log_file


def load_group(h5path: str):
    """Load the merged run's h5 file."""
    print ("------------------------------- Loading pools ----------------------------------------------")
    group = ms.load(h5path, raw=False, filter_variants=False, filter_cells=False)
    print("Pools loaded: ", [s.name for s in group.samples], "\n")
    return group


def remove_doublets_and_add_metadata(group, demux_dir: str, sample_metadata_path: str):
    """
    For each pool:
      - add sample_id and doublet_score
      - remove doublets/unassigned cells
      - add extra row attributes from metadata (per sample_id)
    """
    print ("-------------------- Starting the doublet removal process ----------------------------------")
    sample_metadata = pd.read_csv(sample_metadata_path, sep=";")

    attribute_map = {
        'pool': 'pool',
        'gentli_id': 'gentli_id',
        'ID' : 'Sample',
        'group': 'Group',
        'sample_type': 'Sample Type',
        'sample_tissue': 'tissue',
        'sex': 'SEXatBirth',
        'site': 'Recruitment_Site',
        'age_at_onset': 'AGEATONSET',
        'age_at_death': 'AGEATDEATH',
        'race': 'Racial_Category',
        'ethnicity': 'Ethnicity',
        'disease_duration': 'disease_duration',
        'PMI': 'PMI (h:m)',
        'hem':'hem',
        'hemisphere_ab':'Hemisphere'
    }

    for sample in group:
        pool = sample.name.split('-')[0].lower()
        print(pool, "\n\t",  f"Number of cells: {sample.dna.shape[0]}, Number of positions: {sample.dna.shape[1]}")

        # Add metadata (doublets)
        path = f'{demux_dir}/barcodes/{pool}.csv'
        metadata = pd.read_csv(path, index_col=0)

        unreliable_cells_path = f'{demux_dir}/{pool.capitalize()}/unreliable_cells.txt'
        with open(unreliable_cells_path) as f:
            unreliable_cells = [line.strip() for line in f if line.strip()]

        if (len(unreliable_cells) > 0):
            sample.dna = sample.dna.drop(unreliable_cells, value_type="barcode")
            sample.cnv = sample.cnv.drop(unreliable_cells, value_type="barcode")
        
        # Attach row attrs
        sample.dna.add_row_attr('sample_id', np.array(metadata['sample_id'].to_list()))
        sample.cnv.add_row_attr('sample_id', np.array(metadata['sample_id'].to_list()))
        sample.dna.add_row_attr('doublet_score', np.array(metadata['doublet_pred'].to_list()))
        sample.cnv.add_row_attr('doublet_score', np.array(metadata['doublet_pred'].to_list()))

        # Remove doublets/unassigned
        cells_to_remove = (metadata[metadata['doublet_pred'] == 'High'].index.tolist())
        print("\t", f"No. of cells removed = {len(cells_to_remove)}")
        metadata = metadata.drop(cells_to_remove, axis=0)
        sample.dna = sample.dna.drop(cells_to_remove, value_type="barcode")
        sample.cnv = sample.cnv.drop(cells_to_remove, value_type="barcode")
        print("\t", f"Number of cells : {sample.dna.shape[0]}, Number of positions: {sample.dna.shape[1]}")
        print("\n")


        row_attr_data = {key: [] for key in attribute_map}
        # Populate row data for each sample_id
        for sample_id in sample.dna.row_attrs['sample_id']:
            sample_row = sample_metadata[sample_metadata['Sample'] == sample_id].iloc[0]
            for attr_key, column_name in attribute_map.items():
                row_attr_data[attr_key].append(sample_row[column_name])

        # Add row attributes to the assays
        for attr_key, values in row_attr_data.items():
            data_array = np.array(values)
            sample.dna.add_row_attr(attr_key, data_array)
            sample.cnv.add_row_attr(attr_key, data_array)

def split_pools_into_samples(group):
    """Split each pool by sample_id and regroup."""
    print ("------------------------------ Splitting the pools into samples ----------------------------")
    samples_list = []
    for sample in group:
        samples = sample.split("sample_id")  # List of samples
        samples_list.extend(samples)

    group_new = ms.SampleGroup(samples_list)

    # Set the sample name to sample ID 
    for sample in group_new:
        print(list(set(sample.dna.row_attrs['sample_id'])))
        new_name = np.unique(sample.dna.row_attrs['sample_id'])[0]
        sample.rename(new_name)
    print("Samples extracted: ", [s.name for s in group_new.samples], "\n")
    return group_new


def filter_somatic_variants(group_new, whitelist=None):
    """Filter QC variants per sample, take union + whitelist, subset all samples to final_vars."""
    print ("---------------------------- Filtering good quality variants -------------------------------")

    def custom_filt(sample):
        filt_vars = sample.dna.filter_variants(
            min_dp=10,
            min_gq=30,
            vaf_ref=10,
            vaf_hom=90,
            vaf_het=30,
            min_prct_cells=50,   # genotyped for >50% of cells
            min_mut_prct_cells=(100/sample.dna.shape[0]),  # accommodates variants found in even a single cell
            iterations=10
        )
        return filt_vars
    dna_vars = group_new.apply(custom_filt)

    # Check the number of filtered variants
    print("Number of filtered variants")
    for i in range(len(group_new.samples)):
        sample = group_new.samples[i]
        print("\t", f"{sample.name}: {len(dna_vars[i])}")

    # Union of all variants
    var_union = list(set().union(*dna_vars))
    print(f"Number of all variants selected (union): {len(var_union)}")

    if whitelist is not None:
        final_vars = list(set(var_union).union(whitelist))
    else:
        final_vars = list(var_union)

    print(f"Number of all variants selected (union + white_list): {len(final_vars)}", "\n")

    # Number of variants for each sample (before)
    og_num_vars = [s.dna.shape[1] for s in group_new.samples]
    print(f"sample.dna variants dimension before variant filtering : {list(set(og_num_vars))}")
    print("The number of variants selected across each sample are equal:", all(x == og_num_vars[0] for x in og_num_vars))

    # Subset all samples with the same variants
    for sample in group_new:
        sample.dna = sample.dna[:, final_vars]

    # Number of variants for each sample (after)
    new_num_vars = [s.dna.shape[1] for s in group_new.samples]
    print(f"sample.dna variants dimension after variants selection : {list(set(new_num_vars))}")
    print("The number of variants filtered for each sample are equal:", all(x == new_num_vars[0] for x in new_num_vars), "\n")

    # Summary of DNA assay, per sample
    print("Dimensions after variant selection:")
    for sample in group_new:
        print("\t", f"{sample.name}: {sample.dna.shape}")

    return final_vars


def merge_and_save(group_new, out_path: str):
    """Merge samples into one group and save h5."""
    print (" ------------------------------ Merging samples --------------------------------------------")
    combined = group_new.merge()
    ms.save(combined, out_path)
    print(f"Merged samples h5 file exported to {out_path}")


def main():
    h5path  = '/results/rr/study/hg38s/study232-missionbio_TDP-C/data/merged_assays/hg38/merged.h5'
    demux_results_dir = '/results/rr/study/hg38s/study232-missionbio_TDP-C/data/final/hg38/demultiplexing_results'
    sample_metadata_path = '/results/rr/study/hg38s/study232-missionbio_TDP-C/data/metadata_hemisphere_PMI.csv'
    log_path = '/results/rr/study/hg38s/study232-missionbio_TDP-C/data/outputs/log_files/hg38/multisample_integration_pools_allvars_vafref10_vafhet30_vafhom90.log'
    out_path = '/results/rr/study/hg38s/study232-missionbio_TDP-C/data/outputs/integrated_h5_files/hg38/combined_samples_pools_allvars_vafref10_vafhet30_vafhom90.h5'
    log_file = setup_logging(log_path)

    group = load_group(h5path)
    remove_doublets_and_add_metadata(group, demux_results_dir, sample_metadata_path)
    print("Doublets and unreliable cells removed, metadata added")
    group_new = split_pools_into_samples(group)
    print("Samples successfully split")
    filter_somatic_variants(group_new)
    print("Variant filtering complete")
    merge_and_save(group_new, out_path)
    print("h5 file exported")
    log_file.close()

if __name__ == "__main__":
    main()