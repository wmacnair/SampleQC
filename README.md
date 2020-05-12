# SampleQC

:wave: Hey! :D :wave:

`SampleQC` is an R package for robust multivariate, multi-celltype, multi-sample quality control for single cell RNA-seq. QC is typically done by looking at measures of experimental quality for each cell (such as library size, number of features observed, mitochondrial proportion), and identifying cells which are outliers in some way (e.g. high mitochondrial proportion or small library size).

Standard QC approaches suffer from various flaws:

* Outlier detection assuming that the sample consists of only one celltype introduces biases in the cells excluded. For example, excluding cells with low library sizes will preferentially exclude cell _types_ with small library sizes, even though these cells may be perfectly healthy.
* Univariate outlier detection misses cells which are not outliers with respect to any individual QC metric, but are clearly outliers with respect to the observed distribution.
* Identifying cell outliers within each sample individually means that if a celltype is present in only a small proportion in a given sample, it may be identified as an outlier and excluded, although it is healthy.
* They do not provide an automatic assessment of quality at the _sample_ level.

Our method `SampleQC` addresses these problems by fitting a flexible Gaussian mixture model to user-selected QC metrics. 


# Installation

To use this development version of the package, run the following lines in R:
```R
# preliminaries
install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
BiocManager::install("BiocStyle")

# install this repo
devtools::install_github('wmacnair/SampleQC', auth_token='5800e0926a5fa148d3f712d100cc71f6b1e71ea9')

# set up libraries
library('SampleQC')
```

This should load all of the code.

Currently the package is in alpha phase, i.e. it's definitely not ready yet. But it should still be pretty usable, and it would be great to get feedback on it.


# Basic use

## Overview

`SampleQC` can be broken down into two main parts:

1. Calculate differences between sample distributions, and plot these to identify groupings of _samples_ (e.g. poor quality samples, batch effects)
2. Within each sample grouping, fit a statistical model to identify cells which are outliers, which gives you something like the standard `outlier`/`good` outputs from other QC approaches.

## Setting up

There are two possible startpoints for `SampleQC`: from a `SingleCellExperiment` object _sce_, or a `data.frame` where you've already calculated QC metrics for all your cells, _qc_df_. 

If you already have a data.frame, it should contain the following columns:

* `cell_id`, a unique identifier for each cell
* `sample_id`, identifier for experimental sample
* at least one QC metric (user's choice, e.g. _log_counts_, _log_feats_, _logit_mito_)
* optionally, some additional annotations of the samples, e.g. _batch_id_, _patient_id_, _condition_

The optional annotations are for checking whether there are differences in sample-level statistics which are consistent across some groups, e.g. _log_counts_ is consistently low in a particular batch.

For both of these options, you then run a quick preparation function as follows:
```R
# either: sce option
qc_dt   = make_qc_dt(sce)

# or: data.frame option
qc_dt   = make_qc_dt(qc_df)

# you can also specify your own qc_names, but these need to be present in your qc_df
qc_dt   = make_qc_dt(qc_df, qc_names=c('log_counts', 'log_feats', 'logit_mito'))
```

An important point about the QC metrics you use: `SampleQC` uses a Gaussian distribution to identify cell outliers. This means that you can't use things like _mito_prop_ directly, but have to transform them first, e.g. via the inverse logistic function. `SampleQC` does this transformation for you, into the variable _logit_mito_ in the _qc_dt_ object.

## Running `SampleQC`

We first define some parameters that we want to use.

```R
# which QC metrics do we want to use? (the most important bit)
qc_names        = c('log_counts', 'log_feats', 'logit_mito')
# which discrete-valued variables do we want to annotate the samples with?
annot_discrete  = c('well_id', 'patient_id', 'condition')
# which continuous-valued variables do we want to annotate the samples with?
annot_cont      = NULL
```

`SampleQC` generates some additional automatic annotations, such as median mitochondrial proportion by sample and sample size.

We then calculate distances between all the samples, and embed this matrix via a couple of dimensionality reduction options.

```R
mmd_list    = calculate_sample_to_sample_MMDs(qc_dt, qc_names, subsample=200, n_times=20, n_cores=16)
mmd_list    = embed_sample_to_sample_MMDs(mmd_list, qc_dt, annot_discrete, annot_cont, n_nhbrs=5)
print(table(mmd_list$mmd_clusts))
```

Next we fit a Gaussian mixture model to each of the sample groupings that we found. The user needs to specify how many clusters to fit in each group of samples. The quickest way to do this is to start with `K=1` for each cluster, plot the results, and then inspect the outputs to find which value is best for each cluster. Fitting to real data can be difficult; you may to try multiple different values of K, and maybe also tweak some parameters.

```R
em_list     = fit_sampleQC(mmd_list, qc_dt, qc_names, K_list=c(1,1,1,1))
# em_list     = fit_sampleQC(mmd_list, qc_dt, qc_names, K_list=c(2,3,2,2))
```

Once you've fit everything, you can render an html report and check whether it makes sense. 

```R
# define you project name and where you want to save these reports
proj_name   = 'my_project'
save_dir    = '/home/work/my_project/qc/'
dir.create(save_dir)

# render the report!
make_SampleQC_report(mmd_list, em_list, save_dir, proj_name)
```


# Future development

Things I'm planning to include:

* Unit tests
* Some kind of testing of the *Rcpp* bits

Things I'm considering including:

* Functions for calculating QC metric values e.g. from an `sce` object. This can be slow; it could be worth trying to speed this up, or at least check if someone else has solved this problem already.
* Some way of checking whether the clustering of samples worked properly (specifically, feeding back any signs of bad fit at the GMM step into the earlier stage).


#  Bugs / Suggestions / Thoughts

Please add anything you like to the _Issues_ page. All feedback enthusiastically received.

Cheers

Will




## System requirements

`SampleQC` requires `R >= 4.0.0` (although it may not make much difference if you use `R 3.6.1`.

