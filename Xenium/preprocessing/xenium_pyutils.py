import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib as mpl
import seaborn as sns
import matplotlib.pyplot as plt
#import squidpy as sq
import pickle
import anndata as ad

from PIL import Image
import scipy
#import cv2
#############################
######### Processing#########
#############################

## load xenium from feature matrix, cells and transcripts
# output from xenium onboard analysis
# feature_mtx.h5, cells.parquet, transcripts.parquet
def loading_xenium(cell_feature_mtx, cells, transcripts):

    adata = sc.read_10x_h5(cell_feature_mtx)
    cells = pd.read_parquet(cells)
    cells.set_index(adata.obs_names, inplace=True)
    adata.obs = cells.copy()
    transcripts = pd.read_parquet(transcripts)
    adata.obsm["spatial"] = adata.obs[["x_centroid", "y_centroid"]].copy().to_numpy()

    sc.pp.filter_cells(adata, min_counts=10)
    sc.pp.filter_genes(adata, min_cells=5)

    return adata

## adata without normalization, raw counts
def normalizing(adata):
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata,target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers['lognormal'] = adata.X.copy()
    return adata

## clustering better in seadragon, ok if not many cells
# use default parameters
def clustering(adata):
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)
    # sc.tl.louvain(adata)
    return adata

## method is wilcoxon, can also try t-test
# group is the clustering method, louvain or leiden
# after ranking genes, data should be saved in pickle
def ranking_gene(adata, group):
    sc.tl.rank_genes_groups(
        adata, groupby=group, method="wilcoxon", key_added=group+"_wilcoxon_rank"
    ) ## default to use adata.raw, after QC, before normalization
    sc.tl.filter_rank_genes_groups(
        adata,
        key=group+"_wilcoxon_rank",
        key_added=group+"_wilcoxon_rank_filtered",
    )
    return adata

# input adata without normalization
def subtyping(adata, batch_key, group='louvain'):
    adata_subtype = normalizing(adata)
    adata_subtype = integrating_harmony(adata, batch_key)
    adata_subtype = ranking_gene(adata,group)
    return adata_subtype

#############################
######### Integration########
#############################

## harmony integration
# log normalized h5ad and the batch key
def integrating_harmony(adata, key):
    sc.tl.pca(adata)
    sc.external.pp.harmony_integrate(adata, key = key)
    sc.pp.neighbors(adata, use_rep='X_pca_harmony')
    # sc.tl.louvain(adata)
    sc.tl.leiden(adata)
    sc.tl.umap(adata)
    return adata


#############################
####### Visualization########
#############################

def plotting_rankedgenes(adata, group, out_dir=None):
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby=group,
        standard_scale="var",
        n_genes=6,
        key=group+"_wilcoxon_rank",
        # save=save1,
        show = False)
    if out_dir is not None:
        plt.savefig(f'{out_dir}/wilcoxon_rank_{group}.pdf')
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby=group,
        standard_scale="var",
        n_genes=6,
        key=group+"_wilcoxon_rank_filtered",
        # save=save2,
        show = False)
    if out_dir is not None:
        plt.savefig(f'{out_dir}/wilcoxon_rank_filtered_{group}.pdf')



def plotting_overlay(hne, adata, alpha = 0.7, color='louvain',group=None):
    '''
    hne: path to the hne image,
    level: image layer extracted from xenium morphology_mip.ome
    '''
    fig, ax = plt.subplots()

    image = Image.open(hne)
    scale_factor = 0.85  # Resize to half of the original size

    # Calculate the new width and height
    new_width = int(image.width * scale_factor)
    new_height = int(image.height * scale_factor)

    # Resize the image using the new width and height
    resized_image = image.resize((new_width, new_height))
    ax.imshow(resized_image)
    sq.pl.spatial_scatter(adata,
                          shape = None,
                            color=color,
                            groups = group,
                            wspace = 0.1,
                            alpha = alpha,
                            size = 0.5,
                            ax = ax)

import subprocess
def bsub(
        SRCDIR,
        JOBNAME,
        COMMAND,
        QUEUE='e40medium',
        WALLT='24:00',
        NNODE=12,
        MEMORY=300,

         ):

    job_sub = f"bsub \
        -J {JOBNAME} \
        -o {SRCDIR}/{JOBNAME}.o.txt \
        -cwd {SRCDIR} \
        -q {QUEUE}\
        -W {WALLT} \
        -n {NNODE} \
        -M {MEMORY} \
        -R rusage[mem={MEMORY}] \
        /bin/bash -c '{COMMAND}'"
    subprocess.Popen(job_sub, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE )


