import squidpy as sq
import pandas as pd
import numpy as np
import os, sys
import scanpy as sc
import anndata as ad
import cellcharter as cc
import argparse
import warnings
warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description='Run CellCharter spatial clustering on an h5ad file.'
)
parser.add_argument('--input',       required=True,  type=str,
                    help='Path to input .h5ad file')
parser.add_argument('--output',      required=True,  type=str,
                    help='Path to output .h5ad file')
parser.add_argument('--sample_key',  default='sample', type=str,
                    help='obs column used as sample/library key (default: sample)')
parser.add_argument('--batch_key',   default='sample',   type=str,
                    help='obs column used as batch key for scVI (only used when --use_rep scvi)')
parser.add_argument('--use_rep',     default='X_pca_harmony', type=str,
                    help='Embedding to use. Pass "scvi" to train scVI first; '
                         'otherwise any key in obsm (default: X_pca_harmony)')
parser.add_argument('--scvi_integration', action='store_true',
                    help='Whether to run clustering with scVI embedding directly')
parser.add_argument('--n_layers',    default=3,      type=int,
                    help='Number of neighbor aggregation layers for CellCharter (default: 3)')
parser.add_argument('--n_clusters',  default=None,   type=int,  nargs='+',
                    help='One or more cluster numbers, e.g. --n_clusters 7 8 9 10. '
                         'Each k gets its own obs column. If omitted, auto-K selection is used.')
parser.add_argument('--autok_min',   default=2,      type=int,
                    help='Minimum k for auto-K search (default: 2)')
parser.add_argument('--autok_max',   default=10,     type=int,
                    help='Maximum k for auto-K search (default: 10)')
parser.add_argument('--autok_runs',  default=10,     type=int,
                    help='Number of runs per k for auto-K stability (default: 10)')
parser.add_argument('--gpus',        default=[0],    type=int, nargs='+',
                    help='One or more GPU indices (default: 0). '
                         'Single k + multiple GPUs → DDP across all GPUs. '
                         'Multiple k + multiple GPUs → each k runs on its own GPU in parallel. '
                         'Set to -1 for CPU.')
parser.add_argument('--seed',        default=12345,  type=int,
                    help='Random seed (default: 12345)')
parser.add_argument('--out_key',     default='spatial_cluster', type=str,
                    help='obs column name for the resulting cluster labels (default: spatial_cluster)')
args = parser.parse_args()

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
_use_cpu    = len(args.gpus) == 1 and args.gpus[0] < 0
accelerator = 'cpu' if _use_cpu else 'gpu'


def train_scvi(adata):
    import scvi
    scvi.settings.seed = args.seed
    setup_kwargs = dict(layer='counts')
    if args.batch_key is not None:
        setup_kwargs['batch_key'] = args.batch_key
    scvi.model.SCVI.setup_anndata(adata, **setup_kwargs)
    model = scvi.model.SCVI(adata)
    # scVI DDP requires NCCL which needs a modern driver; use single GPU to be safe
    scvi_devices = 1 if _use_cpu else [args.gpus[0]]
    print(f'Training scVI on device(s) {scvi_devices} ...')
    model.train(
        early_stopping=True,
        plan_kwargs=dict(optimizer='AdamW', reduce_lr_on_plateau=True),
        accelerator=accelerator,
        devices=scvi_devices,
    )
    adata.obsm['X_scVI'] = model.get_latent_representation(adata).astype(np.float32)
    print('scVI training complete. Embedding stored in obsm["X_scVI"].')
    return adata

def scvi_integration(adata, use_rep="X_scVI"):
    # Use key_added/neighbors_key so results are stored under scVI-specific keys
    # and do not overwrite harmony neighbors, leiden, or umap
    n_pcs = adata.obsm[use_rep].shape[1]
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=15, n_pcs=n_pcs,
                    key_added='neighbors_scVI')
    sc.tl.leiden(adata, neighbors_key='neighbors_scVI', key_added='leiden_scVI')
    # Rename existing harmony UMAP before sc.tl.umap overwrites X_umap
    if 'X_umap' in adata.obsm:
        adata.obsm['X_umap_harmony'] = adata.obsm.pop('X_umap')
    sc.tl.umap(adata, neighbors_key='neighbors_scVI')
    adata.obsm['X_umap_scVI'] = adata.obsm.pop('X_umap')
    return adata

def build_spatial_graph(adata, embedding_key):
    sq.gr.spatial_neighbors(
        adata,
        library_key=args.sample_key,
        coord_type='generic',
        delaunay=True,
    )
    cc.gr.remove_long_links(adata)
    cc.gr.aggregate_neighbors(
        adata,
        n_layers=args.n_layers,
        use_rep=embedding_key,
        out_key='X_cellcharter',
        sample_key=args.sample_key,
    )


def _make_trainer_params(gpu_list):
    """Build trainer_params for a given list of GPU indices (or CPU)."""
    if _use_cpu:
        devices = 1
    else:
        devices = gpu_list
    params = dict(accelerator=accelerator, devices=devices)
    if not _use_cpu and len(gpu_list) > 1:
        params['strategy'] = 'ddp'
    return params



def cluster_fixed_k(adata, k, gpu_list):
    out_key = args.out_key if len(args.n_clusters) == 1 else f'{args.out_key}_k{k}'
    trainer_params = _make_trainer_params(gpu_list)
    print(f'Clustering k={k} on GPU(s) {gpu_list}, labels → obs["{out_key}"] ...')
    gmm = cc.tl.Cluster(n_clusters=k, random_state=args.seed, trainer_params=trainer_params)
    gmm.fit(adata, use_rep='X_cellcharter')
    adata.obs[out_key] = gmm.predict(adata, use_rep='X_cellcharter')
    print(f'[k={k}] Done.')


def cluster_autok(adata):
    print(f'Auto-K clustering over k={args.autok_min}–{args.autok_max}, '
          f'{args.autok_runs} runs each ...')
    autok = cc.tl.ClusterAutoK(
        n_clusters=(args.autok_min, args.autok_max),
        max_runs=args.autok_runs,
        model_params=dict(
            random_state=args.seed,
            trainer_params=_make_trainer_params(args.gpus),
        ),
    )
    autok.fit(adata, use_rep='X_cellcharter')
    adata.obs[args.out_key] = autok.predict(adata, use_rep='X_cellcharter')
    # Save stability figure alongside output
    fig_path = os.path.splitext(args.output)[0] + '_autok_stability.png'
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cc.pl.autok_stability(autok)
    plt.savefig(fig_path, bbox_inches='tight', dpi=150)
    plt.close()
    print(f'Auto-K stability plot saved to {fig_path}')
    print(f'Done. Cluster labels stored in obs["{args.out_key}"].')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def _patch_null_reader():
    """Register a null-type reader for anndata versions that lack it.
    Needed when h5ad files store None values with encoding_type='null'."""
    try:
        import h5py
        from anndata._io.specs import IOSpec
        from anndata._io.specs.registry import IORegistryError
        # Try reading a null IOSpec — if it fails, register a handler
        try:
            from anndata._io.specs.registry import _REGISTRY
            _REGISTRY.get_read(h5py.Dataset, IOSpec('null', '0.1.0'))
        except (IORegistryError, Exception):
            from anndata._io.specs import _REGISTRY as REG
            @REG.register_read(h5py.Dataset, IOSpec('null', '0.1.0'))
            def read_null(elem, _reader=None):
                return None
    except Exception:
        pass  # best-effort; will surface original error if it still fails


if __name__ == '__main__':
    _patch_null_reader()
    print(f'Reading {args.input} ...')
    adata = sc.read_h5ad(args.input)

    # Ensure spatial coordinates exist
    if 'spatial' not in adata.obsm:
        if 'centroid_x' in adata.obs.columns and 'centroid_y' in adata.obs.columns:
            adata.obsm['spatial'] = adata.obs[['centroid_x', 'centroid_y']].to_numpy()
        else:
            raise ValueError(
                'No "spatial" key in obsm and no centroid_x/centroid_y columns in obs.'
            )

    # Determine embedding key
    if args.use_rep.lower() == 'scvi':
        if 'counts' not in adata.layers:
            adata.layers['counts'] = adata.X.copy()
        adata = train_scvi(adata)
        if args.scvi_integration:
            adata = scvi_integration(adata)
        embedding_key = 'X_scVI'
    else:
        embedding_key = args.use_rep
        if embedding_key not in adata.obsm:
            raise KeyError(
                f'Embedding "{embedding_key}" not found in obsm. '
                f'Available: {list(adata.obsm.keys())}'
            )

    print(f'Using embedding: {embedding_key}')
    print(f'Building spatial graph (n_layers={args.n_layers}) ...')
    build_spatial_graph(adata, embedding_key)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    if args.n_clusters is not None:
        n_gpus = len(args.gpus)
        multi_k = len(args.n_clusters) > 1
        for i, k in enumerate(args.n_clusters):
            # Round-robin GPU assignment across k values; single k uses all GPUs (DDP)
            gpu_list = [args.gpus[i % n_gpus]] if (multi_k and not _use_cpu) else args.gpus
            cluster_fixed_k(adata, k, gpu_list)
    else:
        cluster_autok(adata)

    print(f'Writing output to {args.output} ...')
    adata.write_h5ad(args.output)
    print('All done.')
