## run banksy for p08R1
## banksy neighborhood
import scanpy as sc
import sys
sys.path.append('/Users/bzhao2/Library/CloudStorage/OneDrive-InsideMDAnderson/Projects-AkdemirLab/pipelines/Banksy_py')
import banksy
from banksy_utils.color_lists import spagcn_color
from banksy.initialize_banksy import initialize_banksy
from banksy.run_banksy import run_banksy_multiparam
from banksy_utils.color_lists import spagcn_color
from banksy.labels import Label
import argparse
import pickle
import os
import numpy as np

# sys.path.append('/Users/bzhao2/Library/CloudStorage/OneDrive-InsideMDAnderson/Projects-AkdemirLab/pipelines/PyUtils')
# import xenium_pyutils

def run_banksy(data,ncluster, lamb=0.8):
    bdata = data.copy()
    # bdata = adata[adata.obs['sample'] == sample]
    resolutions = None # automatically searching matched resolution to max_labels
    pca_dims = [20] # number of dimensions to keep after PCA
    lambda_list = [lamb] # lambda
    k_geom = 15 # 15 spatial neighbours
    max_m = 1 # use AGF
    coord_keys = ('x', 'y', 'spatial') # spatial coordinates in adata.obs
    annotation_key = 'cell_type' # cell type annotation
    cluster_algorithm = 'leiden'
    output_dir = '/rsrch4/scratch/genomic_med/bzhao2/akdemirlab-projects/01_catalyst/scripts/NG_review/Oct212025'

    banksy_dict = initialize_banksy(
    bdata,
    coord_keys,
    k_geom,
    nbr_weight_decay= "scaled_gaussian",
    max_m=max_m,
    plt_edge_hist=False,
    plt_nbr_weights=False,
    plt_agf_angles=False,
    plt_theta=False)

    results_df = run_banksy_multiparam(
            bdata,
            banksy_dict,
            lambda_list,
            resolutions,
            color_list = spagcn_color,
            max_m = max_m,
            filepath = output_dir,
            key = coord_keys,
            pca_dims = pca_dims,
            annotation_key = annotation_key,
            max_labels = ncluster,
            cluster_algorithm = cluster_algorithm,
            match_labels = False, ## can be true to compare different parameters
            savefig = False,
            add_nonspatial = False,
            variance_balance = False,
        )
    bdata.obs[f'banksy_nhood_{ncluster}'] = results_df.loc[results_df.index[0], 'labels'].dense.tolist()
    return bdata.obs[['cell_id', f'banksy_nhood_{ncluster}']].set_index('cell_id').to_dict()[f'banksy_nhood_{ncluster}']

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run banksy with configurable parameters.")
    # parser.add_argument("-i", "--input", required=True, help="Input .h5ad file")
    parser.add_argument("-s", "--sample", default=None, help="If provided, subset by obs['sample']==VALUE before running")
    parser.add_argument("-l", "--lambda", dest="lamb", type=float, default=0.8, help="Lambda regularization value (default: 0.8)")
    # parser.add_argument("-o", "--output", default=None, help="Output .h5ad path (default: input with _banksy appended)")
    args = parser.parse_args()
    # load full dataset, optionally subset for running banksy
    adata_full = sc.read_h5ad('/rsrch4/scratch/genomic_med/bzhao2/akdemirlab-projects/01_catalyst/scripts/NG_review/Oct082025/data/all_adata_combined_excluding_p62R1_harmony_more_filtering_2.h5ad')
    adata_full.obs['x'] = adata_full.obsm['spatial'][:,0]
    adata_full.obs['y'] = adata_full.obsm['spatial'][:,1]
    adata_run = adata_full[adata_full.obs["sample"] == args.sample].copy()

    n_cluster_nhood ={}
    for n_cluster in range(2, 10):
        n_cluster_nhood[f'{n_cluster}'] = run_banksy(adata_run, n_cluster)
    
    # save n_cluster_nhood dict to pickle
    out_dir = '/rsrch4/scratch/genomic_med/bzhao2/akdemirlab-projects/01_catalyst/scripts/NG_review/Oct212025/banksy_results'
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f'{args.sample}_banksy.pkl')
    with open(out_path, 'wb') as f:
        pickle.dump(n_cluster_nhood, f, protocol=pickle.HIGHEST_PROTOCOL)
