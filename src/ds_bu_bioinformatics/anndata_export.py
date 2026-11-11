"""Lossless samples-by-features H5AD export with clinical and feature annotations."""
import json

import anndata as ad
import numpy as np
import pandas as pd


def write_anndata(expression, clinical, genes, metadata, path):
    obs = clinical.set_index('sample_id').loc[expression.columns].copy()
    obs.index.name = 'sample_id'
    var = genes.set_index('feature_id').loc[expression.index].copy()
    # Categoricals preserve missing strings without requiring experimental H5AD
    # nullable-string support. Nullable floats use IEEE NaN for missing values.
    for frame in (obs, var):
        for col in frame:
            if pd.api.types.is_string_dtype(frame[col]) or frame[col].dtype == object:
                frame[col] = pd.Categorical(frame[col].astype(object))
            elif isinstance(frame[col].dtype, pd.Float64Dtype | pd.Float32Dtype):
                frame[col] = frame[col].to_numpy(dtype=float, na_value=np.nan)
    data = ad.AnnData(X=expression.to_numpy().T.copy(), obs=obs, var=var,
                      uns={'metadata_json': json.dumps(metadata, sort_keys=True)})
    data.write_h5ad(path, compression='gzip')
    saved = ad.read_h5ad(path)
    if saved.obs_names.tolist() != expression.columns.tolist() or saved.var_names.tolist() != expression.index.tolist():
        raise ValueError('AnnData export changed sample/feature order')
    if not np.array_equal(saved.X, expression.to_numpy().T):
        raise ValueError('AnnData export changed deposited expression values')
    pd.testing.assert_frame_equal(saved.obs, obs, check_dtype=False, check_categorical=False)
    pd.testing.assert_frame_equal(saved.var, var, check_dtype=False, check_categorical=False)
    if saved.uns['metadata_json'] != data.uns['metadata_json']:
        raise ValueError('AnnData export changed provenance metadata')
