#--------------------------------------------------------------------------------
# filename : harmony_integration.py
# Date : 2024-10-02
# contributor : Bo Zhao
# function: integrate xenium data
# input: merged h5ad
# output: harmony integrated h5ad
#--------------------------------------------------------------------------------




import scanpy as sc
import pandas as pc
import numpy as np
import xenium_pyutils
import random
random.seed(42)
np.random.seed(42)


import argparse
parser = argparse.ArgumentParser(description='Harmony integration')
parser.add_argument('--input_file', type=str, help='Input file')
parser.add_argument('--output_file', type=str, help='Output directory')
parser.add_argument('--n_highly_variable', type=int, default=2000, help='Number of highly variable genes')
parser.add_argument('--batch_key', type=str, default='sample', help='Batch key for harmony integration')
parser.add_argument('--exclude_samples', nargs='+', default=None, help='Batch key for harmony integration')
# parser.add_argument('--xenium', type=bool, default=True, help='Batch key for harmony integration')
args = parser.parse_args()


def run_harmony(input, output, batch, ntop):

    adata = sc.read_h5ad(input)
    if args.exclude_samples != None:
        adata = adata[~adata.obs['sample'].isin(args.exclude_samples)]
    adata = xenium_pyutils.normalizing(adata)

    # if not args.xenium:
    #     sc.pp.highly_variable_genes(adata, n_top_genes=ntop)

    adata = xenium_pyutils.integrating_harmony(adata, batch)
    adata.write_h5ad(output)

if __name__ == "__main__":
    run_harmony(args.input_file, args.output_file, args.batch_key, args.n_highly_variable)

