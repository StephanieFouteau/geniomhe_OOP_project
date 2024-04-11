# scanpy_geniomhe
Python Bioinformatics code for single-cell omics data analysis

## Goal

The goal of this modification of scanpy is to allow some *dynamic plots* in order to have more information, and to generate a *QC report* to give a summary of the results after the different important preprocessing steps.

## Instructions/use

### Installation

To use our library, you have to import our git. Our functions are in the folder `oopsc`.

First, you have to install dependencies, like this :

```bash
git clone git@github:StephanieFouteau/geniomhe_OOP_project.git oopsc
cd oopsc
pip install -e .
```

### Use

Then, you can use our work like anyother library.

```python
import oopsc
``` 

### Dynamics plots

The dynamic plot functions, `dynamic_highly_variable_genes` and `dynamic_scatter`, do not make all plots of this type dynamic, but give plots for specific preprocessing steps:

- one for *quality control* to show low-quality cells to remove.
- the other for the *feature selection* : replace the scanpy plot `highly_variable_genes`

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

