import scanpy as sc
import plotly as py
import plotly.express as px


"""Dynamically plot dispersions versus means for genes. 

adata
    Result of :func:`~scanpy.pp.highly_variable_genes`.
   
"""

def dynamic_plot_highly_variable_genes(adata):
    
    fig = px.scatter(adata.var, x="means", y="dispersions", color="highly_variable", hover_data=['gene_ids'])
    fig.show()
    return fig
