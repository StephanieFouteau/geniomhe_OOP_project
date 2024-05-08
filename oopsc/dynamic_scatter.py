import scanpy as sc
import plotly as py
import plotly.express as px

def dynamic_plot_scatter_genes(adata):
    fig_c_genes = px.scatter(adata.obs, x="total_counts", y="pct_counts_mt", color="n_genes_by_counts")
    # fig_c_genes.show()
    return(fig_c_genes)

def dynamic_plot_scatter_mt(adata):
    fig_c_mt = px.scatter(adata.obs, x="total_counts", y="n_genes_by_counts", color="pct_counts_mt")
    # fig_c_mt.show()
    return(fig_c_mt)

def dynamic_plot_scatter_total(adata):
    fig_c_total = px.scatter(adata.obs, x="pct_counts_mt", y="n_genes_by_counts", color="total_counts")
    # fig_c_total.show()
    return(fig_c_total)
