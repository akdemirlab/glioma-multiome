## calculate the distance 
## or mean distance of 8 nearest neighbors

from scipy.spatial import KDTree
import pandas as pd
import numpy as np

def get_nearest_distance_by_cell_type(data, source_cell_type, target_cell_type, mean_distance=False):
    """
    Calculate the distance between cells of source_cell_type and their nearest target_cell_type neighbor.
    
    Parameters:
    -----------
    adata : AnnData
        AnnData object containing spatial coordinates
    source_cell_type : str
        Cell type for which to find nearest neighbors
    target_cell_type : str
        Cell type to find as nearest neighbors
    
    Returns:
    --------
    pd.DataFrame
        DataFrame containing the distances for each source cell
    """
    results = []
    # Get coordinates of source and target cells
    source_cells = data[data.obs['cell_type'] == source_cell_type] # pericyte
    target_cells = data[data.obs['cell_type'] == target_cell_type]
    
    # Extract coordinates
    source_coords = source_cells.obs[['x', 'y']].values
    target_coords = target_cells.obs[['x', 'y']].values
    
    # Build KDTree for efficient nearest neighbor search
    tree = KDTree(target_coords)
    
    k_neighbors = 2 if target_cell_type == source_cell_type else 1
    if mean_distance:
        if target_cell_type == source_cell_type:
            distances_mean, _ = tree.query(source_coords, k=9)
            # Remove self (distance 0) from mean calculation
            distances = distances_mean[:, 1:].mean(axis=1)
        else:
            distances_mean, _ = tree.query(source_coords, k=8)
            distances = distances_mean.mean(axis=1)
    else:
        distances, _ = tree.query(source_coords, k=k_neighbors)
        if target_cell_type == source_cell_type:
            # Remove self (distance 0), take the second nearest
            distances = distances[:, 1]
    
    # Create result dataframe
    cell_df = pd.DataFrame({
        'cell_id': source_cells.obs.index,
        'distance': distances,
    })
    
    results.append(cell_df)
    return pd.concat(results, ignore_index=True)

def grid_ids(df, x="x", y="y", xstep=100, ystep=100):
    """
    Assign grid cell IDs to points for block-stratified permutations.
    Bin edges are spaced by step width (xstep, ystep).
    """
    xvals, yvals = df[x].to_numpy(), df[y].to_numpy()
    eps = 1e-9
    xmin, xmax = xvals.min() - eps, xvals.max() + eps
    ymin, ymax = yvals.min() - eps, yvals.max() + eps
    xb = np.arange(xmin, xmax + xstep, xstep)
    yb = np.arange(ymin, ymax + ystep, ystep)

    nx = len(xb) - 1
    ny = len(yb) - 1

    ix = np.clip(np.digitize(xvals, xb) - 1, 0, nx - 1)
    iy = np.clip(np.digitize(yvals, yb) - 1, 0, ny - 1)

    return pd.Series(list(zip(ix, iy)), index=df.index, name="grid_id")

def shuffle_labels_preserve_hist(df, label="label", block=None, rng=None):
    """
    Shuffle labels while preserving class frequencies.
    Optimized with vectorized operations and pre-allocated arrays.
    """
    if rng is None:
        rng = np.random.default_rng(42)
    
    labels = df[label].to_numpy()
    
    if block is None:
        # Global shuffle - fastest path
        unique, counts = np.unique(labels, return_counts=True)
        pool = np.repeat(unique, counts)
        rng.shuffle(pool)
        return pool
    else:
        # Block-stratified shuffle
        new = np.empty(len(labels), dtype=labels.dtype)
        grouped = df.groupby(block, sort=False)
        
        for g, idx in grouped.groups.items():
            idx_arr = idx.to_numpy()
            sub = labels[idx_arr]
            u, c = np.unique(sub, return_counts=True)
            pool = np.repeat(u, c)
            rng.shuffle(pool)
            new[idx_arr] = pool
        
        return new
    
def downsample_cells(inadata, source_cell_type, target_cell_type, seed=42):
    if target_cell_type != source_cell_type:
        drop_cells = []
        for sample in inadata.obs['sample'].unique():
            sample_adata = inadata[inadata.obs['sample']==sample].copy()
            target_cells = sample_adata[sample_adata.obs['cell_type']==target_cell_type].obs_names # make sure the cell names are unique
            source_cells = sample_adata[sample_adata.obs['cell_type']==source_cell_type].obs_names
            n_target_cells = len(target_cells)
            n_source_cells = len(source_cells)

            if n_target_cells >= n_source_cells:
                sampled_target_cells = np.random.default_rng(seed).choice(target_cells, n_source_cells, replace=False)
                drop_cells.extend(list(set(target_cells) - set(sampled_target_cells)))
            else:
                sampled_source_cells = np.random.default_rng(seed).choice(source_cells, n_target_cells, replace=False)
                drop_cells.extend(list(set(source_cells) - set(sampled_source_cells)))
        return inadata[~inadata.obs_names.isin(drop_cells)].copy()
    else:
        return inadata.copy()


def compute_distances_disturbance(indata, sample, mean_distance=False, max_distance=200, 
                      source_cell_type='MES-like', target_cell_type='Pericyte', disturbance=None,seed=42, 
                      block_size=50, return_adata=False):
    sample_adata = indata.copy()
    if disturbance=='downsample':
        sample_adata = downsample_cells(indata, source_cell_type, target_cell_type)
    elif disturbance=='stratified_shuffling':
        sample_adata = indata.copy()
        df = sample_adata.obs.reset_index()
        df['grid_id'] = grid_ids(df, x='x', y='y', xstep=block_size, ystep=block_size)
        shuffled_labels = shuffle_labels_preserve_hist(df, label='cell_type', block='grid_id', rng=np.random.default_rng(seed))
        sample_adata.obs['cell_type'] = shuffled_labels
    elif disturbance=='global_shuffling':
        sample_adata.obs['cell_type'] = sample_adata.obs['cell_type'][np.random.permutation(sample_adata.n_obs)].values

    
    sample_dist = get_nearest_distance_by_cell_type(sample_adata, source_cell_type=source_cell_type, target_cell_type=target_cell_type, mean_distance=mean_distance)
    sample_dist['disturbance'] = disturbance if disturbance is not None else 'none'
    sample_dist['sample'] = sample
    sample_dist = sample_dist[sample_dist['distance'] <= max_distance]
    if return_adata:
        return sample_dist, sample_adata
    else:
        return sample_dist


def get_neighborhood_composition(adata, target_celltype, k=8, max_distance=200.0):
    """
    Get detailed neighborhood composition statistics for target cells.
    
    Args:
        adata: AnnData object
        target_celltype: Target cell type to analyze
        k: Number of nearest neighbors
        max_distance: Maximum distance threshold
    
    Returns:
        dict with keys:
            'ratios': DataFrame of cell type ratios per target cell
            'counts': DataFrame of cell type counts per target cell
            'summary': Series with mean ratio of each cell type across all target cells
            'n_neighbors': Series with actual number of neighbors per target cell
    """
    
    # Get target cells
    target_mask = adata.obs['cell_type'] == target_celltype
    target_cells = adata.obs.index[target_mask].tolist()
    target_indices = np.where(target_mask)[0]
    
    if len(target_cells) == 0:
        raise ValueError(f"No cells found with cell_type == '{target_celltype}'")
    
    print(f"Analyzing {len(target_cells)} target cells of type '{target_celltype}'")
    
    # Extract centroids and build KDTree
    centroids = adata.obs[['centroid_x', 'centroid_y']].values
    tree = KDTree(centroids)
    
    # Get all unique cell types
    cell_types = sorted(adata.obs['cell_type'].unique().tolist())
    
    # Initialize output dataframes
    neighbor_counts_df = pd.DataFrame(0, index=target_cells, columns=cell_types)
    neighbor_ratios_df = pd.DataFrame(0.0, index=target_cells, columns=cell_types)
    n_neighbors_series = pd.Series(0, index=target_cells)
    
    # For each target cell, find nearest neighbors and count cell types
    for target_idx, target_cell_id in zip(target_indices, target_cells):
        # Query k+1 neighbors (including self)
        distances, indices = tree.query(centroids[target_idx:target_idx+1], k=k+1)
        distances = distances[0]
        indices = indices[0]
        
        # Filter neighbors: exclude self and those beyond max_distance
        valid_neighbors_idx = [
            idx for i, idx in enumerate(indices[1:], 1)
            if distances[i] <= max_distance
        ]
        
        n_neighbors_series[target_cell_id] = len(valid_neighbors_idx)
        
        if not valid_neighbors_idx:
            # If no valid neighbors, set uniform distribution
            neighbor_ratios_df.loc[target_cell_id, :] = 1.0 / len(cell_types)
            neighbor_counts_df.loc[target_cell_id, :] = 0
            continue
        
        # Get cell IDs and types of valid neighbors
        neighbor_cell_ids = adata.obs.index[valid_neighbors_idx].tolist()
        neighbor_celltypes = adata.obs.loc[neighbor_cell_ids, 'cell_type'].values
        
        # Count each cell type
        unique, counts = np.unique(neighbor_celltypes, return_counts=True)
        for celltype, count in zip(unique, counts):
            neighbor_counts_df.loc[target_cell_id, celltype] = count
            neighbor_ratios_df.loc[target_cell_id, celltype] = count / len(valid_neighbors_idx)
    
    # Calculate summary statistics (mean ratio across all target cells)
    summary_ratios = neighbor_ratios_df.mean()
    
    return {
        'ratios': neighbor_ratios_df,
        'counts': neighbor_counts_df,
        'summary': summary_ratios,
        'n_neighbors': n_neighbors_series
    }
