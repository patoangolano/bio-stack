"""Verifica mediana de UMI, ntc_id e genes com CPM medio > 5 nos controles D/E/F."""
import sys

import h5py
import numpy as np
import scipy.sparse as sp

for path in sys.argv[1:]:
    with h5py.File(path, "r") as f:
        g = f["X"]
        shape = tuple(g.attrs["shape"])
        data = g["data"][:]
        indices = g["indices"][:]
        indptr = g["indptr"][:]
        obs = f["obs"]
        ntc = obs["ntc_id"]
        if isinstance(ntc, h5py.Group):
            cats = ntc["categories"][:]
            codes = ntc["codes"][:]
            n_ntc = len(np.unique(codes))
            per = np.bincount(codes)
        else:
            vals = ntc[:]
            u, per = np.unique(vals, return_counts=True)
            n_ntc = len(u)
    X = sp.csr_matrix((data.astype(np.float64), indices, indptr), shape=shape)
    lib = np.asarray(X.sum(axis=1)).ravel()
    integral = bool(np.all(data == np.round(data)))
    cpm = sp.diags(1e6 / lib) @ X
    mean_cpm = np.asarray(cpm.mean(axis=0)).ravel()
    print(path.split("/")[-1], "cells", shape[0], "ntc_ids", n_ntc,
          "per_ntc min/max", per.min(), per.max(),
          "median_umi", np.median(lib), "nnz/cell", round(len(data) / shape[0], 1),
          "integral", integral, "genes_cpm_gt5", int((mean_cpm > 5).sum()),
          "max_lib", lib.max())
