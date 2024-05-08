# scanpy_geniomhe
Python Bioinformatics code for single-cell omics data analysis using scanpy library
https://scanpy.readthedocs.io/en/stable/

## Goal

The goal of this modification of scanpy is to allow some *dynamic plots* in order to have more information, and to generate a *QC report* to give a summary of the results after the different important preprocessing steps.

## Instructions/use

### Installation

To use our library, you have to import our git. Our functions are in the folder `oopsc`.

First, you have to install dependencies, like this :

```bash
git clone https://github.com/StephanieFouteau/geniomhe_OOP_project/
cd geniomhe_OOP_project/
pip install -e .
```

### Use

Then, you can use our work like anyother library.

```python
import oopsc
```
 
### Data

Scanpy is based on an AnnData object. adata stores a data matrix adata.X, annotation of observations adata.obs and variables adata.var as pd.DataFrame https://scanpy.readthedocs.io/en/stable/usage-principles.html. Let’s start by building a basic AnnData object:

```adata = sc.read_10x_mtx(
    'data/filtered_gene_bc_matrices/hg19/',  # the directory with the `.mtx` file
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                              # write a cache file for faster subsequent reading
```

### Dynamics plots

The dynamic plot functions, `dynamic_highly_variable_genes` and `dynamic_scatter`, give plots for specific preprocessing steps:

*feature selection* : replace the scanpy plotting function `highly_variable_genes`
Plot dynamically dispersions versus means for genes

- `dynamic_highly_variable_genes` 



*quality control* to show low-quality cells to remove.
In the file `dynamic_scatter`, there are three functions, with the same usage :

- `dynamic_plot_scatter_genes`
- `dynamic_plot_scatter_mt`
- `dynamic_plot_scatter_total`

*To make all the scatter type plots dynamic, we have to directly modify the scatter function of scanpy, and more precisely the `_scatter_obs` function.*

### QC report

Our most important project is the addition of functions to generate a QC-report. We have a version directly in html, and another using dash.

The output is an HTML page, with all the plots showing the results of the important preprocessing steps.

In files:

- `QC_report.py`
- `QC_report_dash.py`

### Git utilisation, main command :
A little help for our collaborative work

- `git fetch`: Retrieves changes from a remote repository.
- `commit (with a message)`: Records changes to the repository with a descriptive message.
- `git pull`: Fetches changes from a remote repository and merges them into the current branch.
- `git push`: Sends local commits to a remote repository.
- `git merge`: Integrates changes from one branch into another.

In Visual Studio Code, go to "Source Control"

## authors

- Naïa Périnelle
- Stéphanie FOUTEAU

