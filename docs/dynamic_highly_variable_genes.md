## Dynamic highly variable genes

```python
import oopsc.dynamic_highly_variable_genes as hvg
```

### Plot dynamically dispersions versus means for genes.

*Takes the result of `scanpy.pp.highly_variable_genes` (adata) as parameter.*
`scanpy.pp.highly_variable_genes` expects logarithmized data :

```python
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata) # with all default values for parameters
```
*Replace scanpy plotting function:*

```python
sc.pl.highly_variable_genes(adata)
```


```python
hvg.dynamic_highly_variable_genes(adata) 
```

This dynamic scatter plot displays either dispersion and mean value and the gene_Id for each gene. 
