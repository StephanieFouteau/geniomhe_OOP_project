"""
Propose the Quality Control of sc_RNAseq data following the scverse Tutorial and output a QC_report in html format 
https://scverse-tutorials.readthedocs.io/en/latest/notebooks/basic-scrna-tutorial.html#quality-control
"""
# Imports
import scanpy as sc
import plotly as py
import plotly.express as px
import plotly.io as pio
import seaborn as sns
from matplotlib import pyplot as plt
import is_outlier as is_outlier


def create_qc_report(adata):

    ## Highest expressed genes over all cells.
    # Computes, for each gene, the fraction of counts assigned to that gene within a cell.
    # The 25_top genes with the highest mean fraction over all cells are plotted as boxplots.

    sc.pl.highest_expr_genes(adata, n_top=25, save=True) # saving figure to file figures/highest_expr_genes.pdf

    ## Doublet detection
    sc.external.pp.scrublet(adata) # Predict doublets using Scrublet [Wolock19].
    doublet_nbr = adata.obs[adata.obs.predicted_doublet == True].predicted_doublet.count()
    print(f"Detected doublets: {doublet_nbr}")
    print(f"Percentage of doublets: {(doublet_nbr/adata.n_obs) : .2%}")

    # Plot histogram of doublet scores for observed transcriptomes and simulated doublets.                                                                                                                                                                                                               
    sc.external.pl.scrublet_score_distribution(adata, save=True) # saving figure to file figures/scrublet_score_distribution.pdf

    ## Filtering low quality cells
    # Specific gene populations
    # mitochondrial genes
    adata.var["mt"] = adata.var_names.str.startswith("MT-") # annotate the group of mitochondrial genes as 'mt'
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL")) # annotate the group of ribosomal genes as 'ribo'

    # Calculate the QC metrics and the proportions of counts for specific gene populations.
    sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo"], inplace=True, percent_top=[25], log1p=True
    )


    # Filtering outlier using the is_outlier function
    # Filtering log1p_total_counts, log1p_n_genes_by_counts and pct_counts_in_top_25_genes with a threshold of 5 MADs.
    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", 5)
    |   is_outlier(adata, "log1p_n_genes_by_counts", 5)
    |   is_outlier(adata, "pct_counts_in_top_25_genes", 5)
    )
    adata.obs.outlier.value_counts()

    # Filter pct_counts_Mt with 3 MADs. Additionally, filtering cells with a percentage of mitochondrial counts exceeding 8 %.
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
        adata.obs["pct_counts_mt"] > 8
    )
    adata.obs.mt_outlier.value_counts()

    # filtering the AnnData object based on these two additional columns.
    print(f"Total number of cells: {adata.n_obs}")

    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
    print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")

    
    ## Normalization
    # Normalizing to median total counts
    sc.pp.normalize_total(adata)
    # log1p transform
    sc.pp.log1p(adata)

    # Plot counts distribution before and after Normalization
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    p1 = sns.histplot(adata.obs["total_counts"], bins=50, kde=True, fill=True, ax=axes[0])
    axes[0].set_title("Total counts")
    p2 = sns.histplot(adata.obs["log1p_total_counts"], bins=50, kde=True, fill=False, ax=axes[1])
    axes[1].set_title("Normalized counts")
    plt.show()
    fig.savefig('./figures/Normalized counts.png')


    ## Feature selection
    # sc.pp.highly_variable_genes() annotate highly variable genes, expects logarithmized data
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
    # Plot dispersions or normalized variance versus means for genes.
    sc.pl.dynamic_highly_variable_genes(adata)
    #fig = px.scatter(adata.var, x="means", y="dispersions", color="highly_variable", hover_data=['gene_ids'])
    pio.write_html(fig, file="dispersion versus mean.html", auto_open=False) # save plot as html


    ## Create the HTML report

    html_string = '''
    <html>
        <head>
            <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
            <style>body{ margin:0 100; background:whitesmoke; }</style>
        </head>
        <body>
         <h1>Quality Control</h1>
​
            <!-- *** Section 1 *** --->
            <h2>Section 1: Highest expressed genes</h2>
         <iframe width="1000" height="550" frameborder="0" seamless="seamless" scrolling="no" \
    img src="./figures/highest_expr_genes.pdf"></iframe>


​
         <!-- *** Section 2 *** --->
         <h2>Section 2: Normalization</h2>
         <iframe width="1000" height="550" frameborder="0" seamless="seamless" scrolling="no" \
    img src="./figures/Normalized counts.png"></iframe>

​
    <!-- *** Section 3 *** --->
         <h2>Section 3: Feature selection</h2>
         <iframe width="1000" height="550" frameborder="0" seamless="seamless" scrolling="no" \
    src="dispersion versus mean.html"></iframe>

     </body>
    </html>'''

    report = open('/home/svincent/M1_GENIOMHE/POO2/report.html','w')
    report.write(html_string)
    report.close()

    return report