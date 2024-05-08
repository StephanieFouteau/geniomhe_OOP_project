## Dynamic scatter plots

### remove low-quality reads from the dataset

https://scanpy.readthedocs.io/en/latest/api/preprocessing.html

These covariates might correspond to dying cells (cells with broken membranes whose cytoplasmic mRNA has leaked out and therefore only the mRNA in the mitochondria is still present ⇒ low count depth, few detected genes and a high fraction of mitochondrial reads).

- The number of counts per barcode (count depth)
- The number of genes per barcode
- The fraction of counts from mitochondrial genes per barcode

It is, however, crucial to consider the three QC covariates jointly as otherwise it might lead to misinterpretation of cellular signals : 
- Cells with high fraction of mitochondrial counts might for example be involved in respiratory processes.
- Cells with low or high counts might correspond to quiescent cell populations or cells larger in size.

It is therefore preferred to consider multiple covariates when thresholding decisions on single covariates are made. In general, it is advised to exclude fewer cells and be as permissive as possible to avoid filtering out viable cell populations or small sub-populations

In QC, the first step is to calculate the QC covariates or metric. We compute these using the scanpy function sc.pp.calculate_qc_metrics, which can also calculate the proportions of counts for specific gene populations :


```python
# mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-") # annotate the group of mitochondrial genes as 'mt'
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL")) # annotate the group of ribosomal genes as 'ribo'

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, percent_top=[25], log1p=True) # place calculated metrics in adata’s .obs and .var and computing log1p transformed annotations
```
### dynamic_scatter

```python
import oopsc.dynamic_scatter as dsc
```
# plotting dynamically the three covariates:

- the number of genes expressed in the count matrix
- the total counts per cell
- the percentage of counts in mitochondrial genes

*Three different configurations:*

- `dynamic_plot_scatter_genes`
- `dynamic_plot_scatter_total`
- `dynamic_plot_scatter_mt`

*replace scanpy pltotting function:*

```python
sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
```


```python
dsc.dynamic_plot_scatter_mt(adata)
dsc.dynamic_plot_scatter_genes(adata)
dsc.dynamic_plot_scatter_total(adata)
```




