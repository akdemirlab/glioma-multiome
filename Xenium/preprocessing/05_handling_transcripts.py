"""
05_handling_transcripts.py
Post-proseg transcript handling pipeline:

  Step 1 – Concatenate per-sample proseg outputs into a single AnnData dill
  Step 2 – QC filtering, split WT / Mut, save dill + h5ad
  Step 3 – Extract per-sample transcript metadata from filtered WT / Mut h5ad
             NOTE: Step 3 expects h5ad files with leiden-based low-transcript
             cluster removal already applied (suffix _lowTranscriptLeidenClusterRemoved),
             which is produced by a downstream clustering step between Step 2 and Step 3.
"""

import dill
import pickle
import glob
import os
import subprocess
import random

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy import sparse

random.seed(42)
np.random.seed(42)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

DATA_DIR    = os.path.abspath("../../../data")
RESULTS_DIR = os.path.abspath("../../../results")
SCRIPTS_DIR = os.path.abspath("../../../scripts")

proseg_res_dir  = f'{RESULTS_DIR}/resegment/proseg/prior_expan5/proseg'
concat_h5ad_dir = f'{RESULTS_DIR}/resegment/concat_h5ad/proseg_direct_output'

# filenames written / read across steps
CONCAT_RAW_FNAME    = 'catalyst_proseg_prior_expan5um_raw_name_matched.dill'
WT_DILL_FNAME       = 'catalyst_proseg_prior_expan5um_raw_WT.dill'
MUT_DILL_FNAME      = 'catalyst_proseg_prior_expan5um_raw_Mut.dill'
WT_H5AD_FNAME       = 'catalyst_proseg_prior_expan5um_raw_WT_boundary_feature_removed.h5ad'
MUT_H5AD_FNAME      = 'catalyst_proseg_prior_expan5um_raw_Mut_boundary_feature_removed.h5ad'
# produced after downstream leiden-removal step (input to Step 3)
WT_FILTERED_H5AD    = 'WT_catalyst_proseg_prior_expan5um_raw_noBoundary_lowTranscriptLeidenClusterRemoved.h5ad'
MUT_FILTERED_H5AD   = 'Mut_catalyst_proseg_prior_expan5um_raw_noBoundary_lowTranscriptLeidenClusterRemoved.h5ad'
WT_TRANSCRIPTS_PKL  = 'wt_transcripts_proseg.pkl'
MUT_TRANSCRIPTS_PKL = 'mut_transcripts_proseg.pkl'

# QC thresholds (Step 2)
MIN_TRANSCRIPTS = 20
MAX_TRANSCRIPTS = 800
MIN_GENES       = 4

# samples excluded from analysis
XENIUM_BLACKLIST = ['p16R2', 'p25P', 'p33P', 'p38R2', 'p40P', 'p41P',
                    'p20R2', 'p21R1', 'p32R2', 'p33R1', 'p62R1']


# ---------------------------------------------------------------------------
# Step 1: Concatenate per-sample proseg outputs → single AnnData dill
# ---------------------------------------------------------------------------

def load_to_adata(sdir: str) -> ad.AnnData:
    X   = pd.read_csv(f'{sdir}/expected-counts.csv.gz')
    obs = pd.read_csv(f'{sdir}/cell-metadata.csv.gz')

    adata = ad.AnnData(X[sorted([c for c in X.columns if 'UnassignedCodeword' not in c])])
    obs.set_index(adata.obs_names, inplace=True)
    obs['cell'] = [str(x) for x in obs['cell']]
    adata.obs = obs.copy()
    adata.obs['sample'] = os.path.basename(sdir)
    return adata


def concat_proseg_samples(data_dir: str = proseg_res_dir,
                          out_dir: str  = concat_h5ad_dir,
                          fname: str    = CONCAT_RAW_FNAME) -> ad.AnnData:
    """Read all per-sample proseg result subdirectories and concatenate."""
    os.makedirs(out_dir, exist_ok=True)
    samples = subprocess.run(["ls", data_dir], capture_output=True, text=True).stdout.splitlines()

    all_adata = [load_to_adata(os.path.join(data_dir, ss)) for ss in samples]
    adata = ad.concat(all_adata)

    out_path = os.path.join(out_dir, fname)
    with open(out_path, 'wb') as f:
        dill.dump(adata, f)
    print(f"[Step 1] Saved concatenated AnnData → {out_path}")
    return adata


# ---------------------------------------------------------------------------
# Step 2: QC filtering, split WT / Mut, save dill + h5ad
# ---------------------------------------------------------------------------

def filter_and_split(in_dir: str  = concat_h5ad_dir,
                     out_dir: str = concat_h5ad_dir,
                     in_fname: str = CONCAT_RAW_FNAME):
    """Apply QC thresholds, remove blacklisted samples, split WT / Mut."""
    os.makedirs(out_dir, exist_ok=True)

    with open(os.path.join(in_dir, in_fname), 'rb') as f:
        adata = dill.load(f)

    # convert to sparse for cNMF compatibility
    adata = ad.AnnData(
        X   = sparse.csr_matrix(adata.X),
        obs = adata.obs,
        var = adata.var,
    )

    adata.obs['detected_gene_counts'] = (
        np.sum((adata.X >= 1) * 1, axis=1).reshape(-1).tolist()
    )
    adata.obs['proseg_expected_transcript_counts'] = (
        np.sum(adata.X, axis=1).reshape(-1).tolist()
    )
    # mark all genes as highly variable so cNMF uses all genes
    adata.var['highly_variable'] = True

    filtered_samples = sorted(
        set(adata.obs['sample_2'].unique().tolist()) - set(XENIUM_BLACKLIST)
    )
    adata = adata[adata.obs['sample_2'].isin(filtered_samples)]

    fdata = adata[
        (adata.obs['proseg_expected_transcript_counts'] >= MIN_TRANSCRIPTS) &
        (adata.obs['proseg_expected_transcript_counts'] <= MAX_TRANSCRIPTS) &
        (adata.obs['detected_gene_counts'] >= MIN_GENES)
    ]

    adata_wt  = fdata[fdata.obs['idh'] == 'WT']
    adata_mut = fdata[fdata.obs['idh'] != 'WT']

    # save dill (includes multilayer_boundary if present)
    for obj, fname in [(adata_wt, WT_DILL_FNAME), (adata_mut, MUT_DILL_FNAME)]:
        out_path = os.path.join(out_dir, fname)
        with open(out_path, 'wb') as f:
            dill.dump(obj, f)
        print(f"[Step 2] Saved dill → {out_path}")

    # save h5ad without boundary column (causes write errors)
    for obj, fname in [(adata_wt, WT_H5AD_FNAME), (adata_mut, MUT_H5AD_FNAME)]:
        obj = obj.copy()
        if 'multilayer_boundary' in obj.obs.columns:
            del obj.obs['multilayer_boundary']
        out_path = os.path.join(out_dir, fname)
        obj.write_h5ad(out_path)
        print(f"[Step 2] Saved h5ad  → {out_path}")


# ---------------------------------------------------------------------------
# Step 3: Extract per-sample transcript metadata from filtered WT / Mut h5ad
#         Requires leiden-based low-transcript cluster removal to have been
#         applied upstream (files: *_lowTranscriptLeidenClusterRemoved.h5ad)
# ---------------------------------------------------------------------------

def extract_proseg_transcripts(h5ad_dir: str    = concat_h5ad_dir,
                                proseg_dir: str  = proseg_res_dir,
                                out_dir: str     = concat_h5ad_dir):
    """Collect per-sample transcript-metadata.csv.gz, filter to assigned cells,
    and save per-IDH-group dictionaries as pickle files."""
    os.makedirs(out_dir, exist_ok=True)

    wt_adata  = sc.read_h5ad(os.path.join(h5ad_dir, WT_FILTERED_H5AD))
    mut_adata = sc.read_h5ad(os.path.join(h5ad_dir, MUT_FILTERED_H5AD))

    res = subprocess.Popen(
        f'find {proseg_dir} -name transcript-metadata.csv.gz',
        shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8'
    )
    transcripts_csv = [x for x in res.communicate()[0].split('\n') if len(x) > 0]

    # x/y/z are repositioned by proseg; observed_x/y/z are input locations
    wt_transcript_dict  = {}
    mut_transcript_dict = {}

    for ss in transcripts_csv:
        sample  = os.path.basename(os.path.dirname(ss))
        tmp_df  = pd.read_csv(ss)
        tmp_df  = tmp_df[tmp_df['assignment'] != 4294967295]   # unassigned sentinel

        if sample in wt_adata.obs['sample'].cat.categories:
            tmp_adata = wt_adata[wt_adata.obs['sample'] == sample]
            sample_2  = tmp_adata.obs['sample_2'].cat.categories.values[0]
            wt_transcript_dict[sample_2] = tmp_df[
                tmp_df['assignment'].astype('str').isin(tmp_adata.obs['cell'])
            ]
        elif sample in mut_adata.obs['sample'].cat.categories:
            tmp_adata = mut_adata[mut_adata.obs['sample'] == sample]
            sample_2  = tmp_adata.obs['sample_2'].cat.categories.values[0]
            mut_transcript_dict[sample_2] = tmp_df[
                tmp_df['assignment'].astype('str').isin(tmp_adata.obs['cell'])
            ]

    for obj, fname in [
        (wt_transcript_dict,  WT_TRANSCRIPTS_PKL),
        (mut_transcript_dict, MUT_TRANSCRIPTS_PKL),
    ]:
        out_path = os.path.join(out_dir, fname)
        with open(out_path, 'wb') as f:
            pickle.dump(obj, f)
        print(f"[Step 3] Saved transcript dict → {out_path}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    concat_proseg_samples()
    filter_and_split()
    # Step 3: run only after leiden-based low-transcript cluster removal
    # extract_proseg_transcripts()
