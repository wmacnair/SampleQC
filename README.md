# SampleQC

:wave: Hey! :D :wave:

`SampleQC` is an R package for robust multivariate, multi-celltype, multi-sample quality control for single cell RNA-seq. QC is typically done by looking at univariate measures of experimental quality for each cell (such as library size, number of features observed, mitochondrial proportion), and identifying cells that are outliers in some way (e.g. high mitochondrial proportion or small library size).

The main flaw in standard QC approaches that `SampleQC` seeks to correct is that unimodal outlier detection implicitly assumes that the sample consists of only one celltype, introducing biases in the cells excluded. It can preferentially exclude whole cell types with extreme QC metrics (e.g. small library sizes; naturally high mitochondrial proportions), even though these cells may be perfectly healthy. Related to this, identifying cell outliers within each sample individually means that if a celltype is present in only a small proportion in a given sample, it may be identified as an outlier and excluded, although it is healthy.

Our method `SampleQC` addresses these problems by robustly fitting a multivariate Gaussian mixture model, across multiple samples simultaneously. It is intended for large, complex datasets with many samples, although should also work on smaller datasets. By default it uses log counts, log feature-counts and mitochondrial proportions, although users can select which QC metrics they wish to use (including e.g. CITE-seq QC metrics).


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
devtools::install_github('wmacnair/SampleQC')

# set up libraries
library('SampleQC')
```

# Basic use

## Overview

`SampleQC` can be broken down into two main parts:

1. `SampleQC` first calculates dissimilarities between sample distributions, and uses these to identify 'sample groups', i.e. groups of samples with similar distributions. This also provides embeddings of the samples, allowing users to check for batch effects.
2. Within each sample group, `SampleQC` fits a statistical model to identify cells that are potential outliers, based on how unlikely they are under the statistical model.

## Setting up

`SampleQC` starts from a `data.frame` containing the QC metrics for all your cells, _qc_df_. This `data.frame` is easily obtainable from whichever object contains your single cell data. (We also considered writing functions that take single cell container objects as inputs, such as `SingleCellExperiment` or `Seurat` objects. However, using a `data.frame` is more flexibile, and future-proofs `SampleQC` against new types of data and changes to other packages.) 

This `data.frame` should contain the following columns:

* `cell_id`, a unique identifier for each cell
* `sample_id`, identifier for experimental sample
* at least one QC metric (user's choice, e.g. _log_counts_, _log_feats_, and _logit_mito_)
* optionally, some additional annotations of the samples, e.g. _batch_id_, _patient_id_, _condition_

The optional annotations are for checking whether there are differences in sample-level statistics that are consistent across some groups, e.g. _log_counts_ is consistently low in a particular batch.

You then run a quick preparation function. For `SingleCellExperiment` objects, it looks like this:
```R
# calculate QC metrics with `scater`
sce     = scater::addPerCellQC(sce, 
  subsets = list(mito = grep("MT-", rownames(sce))))

# call prep function
qc_dt   = make_qc_dt(colData(sce), 
  sample_var = 'sample_id', 
  qc_names  = c('log_counts', 'log_feats', 'logit_mito'),
  annot_vars  = 'condition'
  )
```

For `Seurat` objects, it looks like this:
```R
# add mito QC metrics
seu[["percent.mt"]] = PercentageFeatureSet(seu, pattern = "^MT-")

# call prep function
qc_dt   = make_qc_dt(seu@meta.data, 
  sample_var  = 'sample_id', 
  qc_names    = c('log_counts', 'log_feats', 'logit_mito'),
  annot_vars  = 'condition'
  )
```

An important point about the QC metrics you use: `SampleQC` uses a mixture of Gaussian distributions to identify cell outliers. This means that the model will fit better for metrics that are more approximately normal. Therefore, metrics like proportions (_percent.mt_) should be transformed first, e.g. via the inverse logistic function. For the variable _percent.mt_, `SampleQC` does this transformation automatically, into the variable _logit_mito_ in the _qc_dt_ object.

## Running `SampleQC`

We first define some parameters that we want to use.

```R
# which QC metrics do we want to use?
qc_names    = c('log_counts', 'log_feats', 'logit_mito')

# which discrete-valued variables do we want to annotate the samples with?
annots_disc = c('well_id', 'patient_id', 'condition')

# which continuous-valued variables do we want to annotate the samples with?
annots_cont = NULL
```

`SampleQC` generates some additional automatic annotations, such as median mitochondrial proportion by sample, sample size, and some others.

We then calculate distances between all the samples, and embed this matrix via dimensionality reduction options. `SampleQC` stores everything neatly in a `SingleCellExperiment` object (mainly in the `colData` entry).

```R
qc_obj    = calc_pairwise_mmds(qc_dt, qc_names, 
  annots_disc = annots_disc, annots_cont = annots_cont, n_cores = 4)
table(colData(qc_obj)$group_id)
```

Next we fit Gaussian mixture models, either one to each of the sample groupings that we found, or to the whole dataset. The user needs to specify how many clusters to fit in each group of samples. The quickest way to do this is to start with `K=1` for each cluster, plot the results, and then inspect the outputs to find which value is best for each cluster. Fitting to real data can require a couple of user iterations, for example trying multiple different values of K.

To fit to each of the sample groupings individually, you use the parameter `K_list`. We recommend first specifying 1 component for each group and rendering a report: `K=1` is extremely quick to fit, and the diagnostic plots in the rendering allow you to check the appropriate number of components for each sample group.

```R
qc_obj    = fit_sampleqc(qc_obj, K_list=rep(1, get_n_groups(qc_obj)))
```

Once the model has fit, you can render an html report and check whether it makes sense. 

```R
# define you project name and where you want to save these reports
proj_name = 'my_project'
save_dir  = '/home/work/my_project/qc/'
dir.create(save_dir)

# render the report
make_sampleqc_report(qc_obj, save_dir, proj_name)
```

This allows you to check whether the number of components for each group looks correct. If not, you can rerun with a different specification of `K_list`:
```R
qc_obj    = fit_sampleqc(qc_obj, K_list=c(2,3,2,2))
```

To fit one model to the whole of the dataset, you use the parameter `K_all`. You might want to do this when you have relatively few samples (e.g. 30 or fewer) and / or when your samples have very similar distributions of QC metrics.

```R
qc_obj    = fit_sampleqc(qc_obj, K_all=2)
```

Once you're happy that the model is identifying outliers correctly, you can extract the outliers found by `SampleQC`by the following:

```R
outliers_dt = get_outliers(qc_obj)
```

There are a couple of ways of tweaking exactly how you get your outliers using the `get_outliers` function. If the plots in the section _Checking for outliers and QC batches via MMD_ of the `SampleQC` report indicate that one group of samples is universally of bad quality, then you can choose to exclude all of those samples via the argument `exc_groups`. If there are groups of cells with sufficient numbers to be modelled as a group by `SampleQC`, but you want to exclude them (for example, a large cluster of cells with extremely high mitochondrial proportions), then you can specify clusters to be excluded with the argument `exc_clusters`. You can see more details by calling `?get_outliers`.

# Preprint

To replicate the analyses, in the preprint, please go to the companion GitHub repo [here](https://github.com/wmacnair/SampleQC_paper_analyses).

#  Bugs / Suggestions / Thoughts

Please add anything you like to the _Issues_ page. All feedback enthusiastically received!

Cheers

Will




## System requirements

`SampleQC` requires `R >= 4.0.0`.

