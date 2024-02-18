from anndata import AnnData
import numpy as np
from scipy.stats import median_abs_deviation
import pandas as pd


def is_outlier(adata: AnnData, metric: str, nmads: int):
    M = adata.obs[metric]
    median = np.median(M)
    mad = median_abs_deviation(M)
    return (M < median - nmads * mad) | (median + nmads * mad < M)


def adaptive_outliers(adata: AnnData, metrics: list[str], nmads: int):
    combined_outliers = pd.Series(False, index=adata.obs.index)

    for metric in metrics:
        current_outliers = is_outlier(adata, metric, nmads)
        combined_outliers = combined_outliers | current_outliers

    return combined_outliers
