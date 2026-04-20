#!/usr/bin/env python3
"""
Efficient parallel processing of spatial transcriptomics regions and clusters.
Outputs convex hull geometry for each cluster identified.
Optimized for processing ~1 million cells.
"""

import argparse
import pickle
import warnings
from pathlib import Path
from multiprocessing import Pool, cpu_count
from functools import partial

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.spatial import KDTree
from scipy.sparse import coo_matrix, csgraph
from sklearn.cluster import DBSCAN
from shapely.geometry import MultiPoint
from tqdm import tqdm
import itertools

warnings.filterwarnings("ignore", category=UserWarning)

# Default constants (can be overridden via command line)
DEFAULT_RADIUS_UM = 1.0
DEFAULT_MIN_POINTS = 3
DEFAULT_MIN_TRANSCRIPTS = 4

# ==================== Core Analysis Functions ====================

def kdtree_clustering(df, radius_um, min_points):
    """
    Find strict spatial regions using KDTree with complete-link or 3-clique criteria.
    Optimized for speed with early exits. Computes convex hull for each cluster.
    
    Args:
        df: DataFrame with x_location, y_location columns
        radius_um: Distance threshold in microns (same as eps_um in DBSCAN)
        min_points: Minimum points per region (same as min_samples in DBSCAN)
    """
    pts = df[["x_location", "y_location"]].to_numpy()
    n = len(pts)
    
    if n < min_points:
        df = df.copy()
        df["cluster_id"] = -1
        return df, [], pd.DataFrame()
    
    kdt = KDTree(pts)
    neigh = kdt.query_ball_tree(kdt, r=radius_um)
    
    # Build adjacency matrix efficiently (upper triangle only)
    rows, cols = [], []
    for i, nb in enumerate(neigh):
        for j in nb:
            if i < j:
                rows.extend([i, j])
                cols.extend([j, i])
    
    if not rows:
        df = df.copy()
        df["cluster_id"] = -1
        return df, [], pd.DataFrame()
    
    A = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n)).tocsr()
    n_comp, labels = csgraph.connected_components(A, directed=False)
    
    df = df.copy()
    df["cluster_id"] = -1
    keep_regions = []
    summaries = []
    
    for comp_id in range(n_comp):
        idx = np.where(labels == comp_id)[0]
        if len(idx) < min_points:
            continue
        
        sub = pts[idx]
        
        # Complete-link criterion
        if len(idx) > 1:
            d = np.sqrt(((sub[:, None, :] - sub[None, :, :]) ** 2).sum(axis=2))
            maxd = d[np.triu_indices_from(d, k=1)].max()
        else:
            maxd = 0.0
        
        is_valid = maxd <= radius_um + 1e-9
        
        # 3-clique criterion if complete-link fails
        if not is_valid and len(idx) >= 3:
            neigh_sets = [set(neigh[i]) for i in idx]
            for i, j, k in itertools.combinations(range(len(idx)), 3):
                gi, gj, gk = idx[i], idx[j], idx[k]
                if (gj in neigh_sets[i] and gk in neigh_sets[i] and
                    gi in neigh_sets[j] and gk in neigh_sets[j] and
                    gi in neigh_sets[k] and gj in neigh_sets[k]):
                    is_valid = True
                    break
        
        if is_valid:
            rid = len(keep_regions)
            df.loc[df.index[idx], "cluster_id"] = rid
            keep_regions.append(idx)
            
            centroid = sub.mean(axis=0)
            
            # Compute convex hull
            try:
                if len(idx) >= 3:
                    hull = MultiPoint(sub).convex_hull
                    hull_wkt = hull.wkt
                else:
                    hull_wkt = None
            except Exception:
                hull_wkt = None
            
            summaries.append({
                "cluster_id": rid,
                "n_transcripts": int(len(idx)),
                "centroid_x_um": float(centroid[0]),
                "centroid_y_um": float(centroid[1]),
                "convex_hull_wkt": hull_wkt,
            })
    
    summary_df = pd.DataFrame(summaries) if summaries else pd.DataFrame()
    return df, keep_regions, summary_df


def dbscan_clustering(df, radius_um, min_points):
    """
    DBSCAN clustering with efficient numpy operations.
    Computes convex hull for each cluster.
    
    Args:
        df: DataFrame with x_location, y_location columns
        radius_um: Distance threshold in microns (eps parameter for DBSCAN)
        min_points: Minimum points per cluster (min_samples parameter for DBSCAN)
    """
    X = df[["x_location", "y_location"]].to_numpy()
    
    if len(X) < min_points:
        df = df.copy()
        df["cluster_id"] = -1
        return df, pd.DataFrame()
    
    labels = DBSCAN(eps=radius_um, min_samples=min_points, metric="euclidean").fit_predict(X)
    df = df.copy()
    df["cluster_id"] = labels
    
    cluster_ids = [c for c in set(labels) if c != -1]
    if not cluster_ids:
        return df, pd.DataFrame()
    
    summaries = []
    for cid in cluster_ids:
        mask = labels == cid
        pts = X[mask]
        centroid = pts.mean(axis=0)
        
        # Compute convex hull
        try:
            if len(pts) >= 3:
                hull = MultiPoint(pts).convex_hull
                hull_wkt = hull.wkt
            else:
                hull_wkt = None
        except Exception:
            hull_wkt = None
        
        summaries.append({
            "cluster_id": cid,
            "n_transcripts": int(len(pts)),
            "centroid_x_um": float(centroid[0]),
            "centroid_y_um": float(centroid[1]),
            "convex_hull_wkt": hull_wkt,
        })
    
    return df, pd.DataFrame(summaries)


# ==================== Parallel Processing ====================

def _process_cell_batch(args):
    """Process a single cell - worker function for parallel execution."""
    cell, cell_data, radius_um, min_points = args
    
    # KDTree regions
    _, _, kdtree_summary = kdtree_clustering(cell_data, radius_um, min_points)
    if not kdtree_summary.empty:
        kdtree_summary = kdtree_summary.copy()
        kdtree_summary['cell_id'] = cell
    else:
        kdtree_summary = pd.DataFrame({'cell_id': [cell]})
    
    # DBSCAN clusters (using same parameters)
    _, dbscan_summary = dbscan_clustering(cell_data, radius_um, min_points)
    if not dbscan_summary.empty:
        dbscan_summary = dbscan_summary.copy()
        dbscan_summary['cell_id'] = cell
    else:
        dbscan_summary = pd.DataFrame({'cell_id': [cell]})
    
    return kdtree_summary, dbscan_summary


def get_region_cluster_summary(sample, cells, transcripts_dict, 
                                min_transcripts=DEFAULT_MIN_TRANSCRIPTS,
                                radius_um=DEFAULT_RADIUS_UM,
                                min_points=DEFAULT_MIN_POINTS,
                                gene='EGFR', 
                                n_processes=None,
                                batch_size=1000):
    """
    Parallel processing of regions and clusters for ~1M cells.
    Outputs convex hull geometry for each cluster.
    
    Args:
        sample: Sample identifier
        cells: List of cell IDs to process
        transcripts_dict: Dictionary of transcript dataframes
        min_transcripts: Minimum EGFR transcripts per cell to include
        radius_um: Distance threshold in microns (used for both KDTree and DBSCAN)
        min_points: Minimum points per region/cluster (used for both methods)
        n_processes: Number of processes (None=auto, uses n_cores-1)
        batch_size: Cells per batch for progress updates
        
    Returns:
        kdtree_summaries_df, dbscan_summaries_df (both with convex_hull_wkt column)
    """
    print(f"Processing sample: {sample}")
    print(f"Parameters: radius={radius_um}μm, min_points={min_points}, min_transcripts={min_transcripts}")
    print(f"Total cells to check: {len(cells):,}")
    
    # Pre-filter EGFR transcripts
    df = transcripts_dict[sample]
    df_egfr = df[df['feature_name'] == gene].copy()
    
    # Filter valid cells (>= min_transcripts)
    cell_counts = df_egfr['cell_id'].value_counts()
    valid_cells = set(cell_counts[cell_counts >= min_transcripts].index)
    cells_to_process = [c for c in cells if c in valid_cells]
    
    print(f"Valid cells (>={min_transcripts} EGFR transcripts): {len(cells_to_process):,}")
    
    if not cells_to_process:
        print("No valid cells to process!")
        return pd.DataFrame(), pd.DataFrame()
    
    # Filter dataframe and group by cell (efficient!)
    df_egfr = df_egfr[df_egfr['cell_id'].isin(cells_to_process)]
    grouped = df_egfr.groupby('cell_id')
    
    # Prepare arguments for parallel processing
    args_list = [
        (cell, group.reset_index(drop=True), radius_um, min_points)
        for cell, group in grouped
    ]
    
    # Determine optimal number of processes
    if n_processes is None:
        n_processes = max(1, cpu_count() - 1)
    
    print(f"Using {n_processes} processes for parallel computation")
    
    # Process with parallel pool
    if n_processes > 1:
        with Pool(n_processes) as pool:
            results = list(tqdm(
                pool.imap(_process_cell_batch, args_list, chunksize=max(1, len(args_list) // (n_processes * 4))),
                total=len(args_list),
                desc="Processing cells",
                unit="cells"
            ))
    else:
        results = [_process_cell_batch(args) for args in tqdm(
            args_list, desc="Processing cells", unit="cells"
        )]
    
    print("Combining results...")
    
    # Separate and concatenate results
    kdtree_summaries = [r[0] for r in results if not r[0].empty]
    dbscan_summaries = [r[1] for r in results if not r[1].empty]
    
    kdtree_df = pd.concat(kdtree_summaries, ignore_index=True) if kdtree_summaries else pd.DataFrame()
    dbscan_df = pd.concat(dbscan_summaries, ignore_index=True) if dbscan_summaries else pd.DataFrame()
    
    print(f"KDTree regions found: {len(kdtree_df):,}")
    print(f"DBSCAN clusters found: {len(dbscan_df):,}")
    
    return kdtree_df, dbscan_df


# ==================== Cell Selection Helper ====================

def get_cells_for_sample(sample, adata, p46p, p67p, selected_celltype=None):
    """
    Extract cell IDs for a given sample and cell types.
    Handles both main adata and special cases for p46P/p67P.
    """
    if sample in ['p46P', 'p67P']:
        adata_sample = p46p if sample == 'p46P' else p67p
        if selected_celltype is not None:
            mask = adata_sample.obs['cell_type_cnmf_transfered'].isin(selected_celltype)
        else:
            mask = np.ones(len(adata_sample), dtype=bool)
        cells = adata_sample[mask].obs['cell_id'].str.split('_').str[1].tolist()
    else:
        if selected_celltype is not None:
            mask = (adata.obs['sample'] == sample) & (adata.obs['cell_type'].isin(selected_celltype))
        else:
            mask = adata.obs['sample'] == sample
        cells = adata[mask].obs['cell_id'].str.split('_').str[1].tolist()
    return cells

# ==================== Main Execution ====================

def main():
    parser = argparse.ArgumentParser(
        description='Process spatial transcriptomics regions and clusters with parallel processing',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--sample', type=str, required=True, 
                        help='Sample name to process')
    parser.add_argument('--data_dir', type=str, default='/path/to/data',
                        help='Directory containing input data')
    parser.add_argument('--res_dir', type=str, default='/path/to/results',
                        help='Directory for output results')
    parser.add_argument('--n_processes', type=int, default=None,
                        help='Number of parallel processes (default: CPU count - 1)')
    parser.add_argument('--min_transcripts', type=int, default=DEFAULT_MIN_TRANSCRIPTS,
                        help='Minimum EGFR transcripts per cell to include in analysis')
    parser.add_argument('--radius_um', type=float, default=DEFAULT_RADIUS_UM,
                        help='Distance threshold in microns (used for both KDTree radius and DBSCAN eps)')
    parser.add_argument('--min_points', type=int, default=DEFAULT_MIN_POINTS,
                        help='Minimum points per region/cluster (used for both KDTree min_size and DBSCAN min_samples)')
    # add gene as an argument if needed in future
    parser.add_argument('--gene', type=str, default='EGFR',
                        help='Gene name to filter transcripts (default: EGFR)')
    parser.add_argument('--selected_types', type=str, nargs='+', default=['AC-like', 'MES-like'],
                        help='List of cell types to include (default: AC-like MES-like)')


    args = parser.parse_args()
    
    # Setup paths
    data_dir = Path(args.data_dir)
    res_dir = Path(args.res_dir)
    res_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*60)
    print("Loading data...")
    print("="*60)
    
    # Load data
    adata = sc.read_h5ad(data_dir / 'new_xenium_samples_filtered_annotated.h5ad')
    p46p = sc.read_h5ad(data_dir / 'p46P_filtered_annotated.h5ad')
    p67p = sc.read_h5ad(data_dir / 'p67P_filtered_annotated.h5ad')
    
    with open(data_dir / 'transcripts_filtered.pkl', 'rb') as f:
        transcripts_dict = pickle.load(f)
    
    # Configuration
    selected_celltype = args.selected_types
    
    # Get cells for sample
    cells = get_cells_for_sample(args.sample, adata, p46p, p67p, selected_celltype)
    
    print(f"\nProcessing sample: {args.sample}")
    print(f"Cell types: {selected_celltype}")
    print(f"Total cells: {len(cells):,}")
    
    print("\n" + "="*60)
    print("Starting parallel analysis...")
    print("="*60 + "\n")
    
    # Process regions and clusters
    kdtree_summary, dbscan_summary = get_region_cluster_summary(
        sample=args.sample,
        cells=cells,
        transcripts_dict=transcripts_dict,
        min_transcripts=args.min_transcripts,
        radius_um=args.radius_um,
        min_points=args.min_points,
        n_processes=args.n_processes, 
        gene=args.gene,
    )
    
    # Save results
    print("\nSaving results...")
    kdtree_file = res_dir / f'{args.sample}_{args.gene}_kdtree_summary.csv'
    dbscan_file = res_dir / f'{args.sample}_{args.gene}_dbscan_summary.csv'
    
    kdtree_summary.to_csv(kdtree_file, index=False)
    dbscan_summary.to_csv(dbscan_file, index=False)
    
    print(f"✓ KDTree results saved: {kdtree_file}")
    print(f"  Columns: cluster_id, n_transcripts, centroid_x_um, centroid_y_um, convex_hull_wkt, cell_id")
    print(f"✓ DBSCAN results saved: {dbscan_file}")
    print(f"  Columns: cluster_id, n_transcripts, centroid_x_um, centroid_y_um, convex_hull_wkt, cell_id")
    print("\n" + "="*60)
    print("Processing complete!")
    print("="*60)


if __name__ == '__main__':
    main()