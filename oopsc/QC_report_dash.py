"""
Propose the Quality Control of sc_RNAseq data following the scverse Tutorial and output a QC_report in html format 
https://scverse-tutorials.readthedocs.io/en/latest/notebooks/basic-scrna-tutorial.html#quality-control
"""
# Imports
from dash import Dash, dcc, html
import scanpy as sc
import plotly as py
import plotly.express as px
import plotly.io as pio
import seaborn as sns
from matplotlib import pyplot as plt
import is_outlier as is_outlier

# The fake data
def prepare_adata():
    results_file = 'write/pbmc3k.h5ad'  # the file that will store the analysis results

    adata = sc.read_10x_mtx(
    'data/filtered_gene_bc_matrices/hg19/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                              # write a cache file for faster subsequent reading

    adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
    return adata

adata = prepare_adata()

# The figures
fig1 = sc.pl.highest_expr_genes(adata, n_top=25, save=True)


# The app
app = Dash(__name__)

app.layout = html.Div(children=[
    html.H1(
        children='QC Report'
    ),

    html.Div(children='Highest expressed genes over all cells.\n\n Computes, for each gene, the fraction of counts assigned to that gene within a cell. The 25_top genes with the highest mean fraction over all cells are plotted as boxplots.'
    ),

    dcc.Graph(
        id='example-graph-2',
        figure=fig1
    )
])

if __name__ == '__main__':
    app.run(debug=True)