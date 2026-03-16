'''
Helper functions for analysis with AnnData objects
'''

import numpy as np
import scipy as sp
import pandas as pd
import anndata as ad
from anndata import AnnData
import h5py
import hdf5plugin


# Function to average gene expression across a group
def grouped_obs_mean(adata, group_key, layer=None, gene_symbols=None):
    if layer is not None:
        getX = lambda x: x.layers[layer]
    else:
        getX = lambda x: x.X
    if gene_symbols is not None:
        new_idx = adata.var[idx]
    else:
        new_idx = adata.var_names

    grouped = adata.obs.groupby(group_key)
    out = pd.DataFrame(np.zeros((adata.shape[1], len(grouped)), dtype=np.float64),columns=list(grouped.groups.keys()),index=adata.var_names)

    for group, idx in grouped.indices.items():
        X = getX(adata[idx])
        out[group] = np.ravel(X.mean(axis=0, dtype=np.float64)).tolist()
    return out