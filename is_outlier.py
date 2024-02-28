"""Single-Cell Analysis in Python.
https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html
"""
# def is-outlier sc best practices
# function that takes a metric, i.e. a column in .obs and the number of MADs (nmad) that is still permissive within the filtering strategy
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier
