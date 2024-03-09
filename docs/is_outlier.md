## Filtering low quality reads


In QC, the first step is to calculate the QC covariates or metric. We compute these using the scanpy function `sc.pp.calculate_qc_metrics`
considering automatic thresholding via MAD (median absolute deviations). The MAD is given by MAD = median(|Xi -median(X)|) with Xi being the respective QC metric of an observation and describes a robust statistic of the variability of the metric. Similar to [Germain et al., 2020](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02136-7), we mark cells as outliers if they differ by 5 MADs which is a relatively permissive filtering strategy. We want to highlight that it might be reasonable to re-assess the filtering after annotation of cells.

we define a function that takes a metric, i.e. a column in .obs and the number of MADs (nmad) that is still permissive within the filtering strategy.